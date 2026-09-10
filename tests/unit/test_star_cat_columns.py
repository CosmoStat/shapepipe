"""The star catalogue's 16 columns are defined twice, and must not drift.

Two writers emit a full_starcat, for two consumers that have to agree about it:

  * ``MergeStarCatPSFEX`` (``src/shapepipe/modules/merge_starcat_package``),
    which the ``merge_starcat`` MODULE RUNNER calls, writing the flat FITS table
    sp_validation opens today;
  * ``workflow/scripts/merge_star_cat.py``, the Snakemake workflow's
    ``star_cat_merge`` rule, writing the per-exposure hdf5 that replaces it
    (CosmoStat/sp_validation#340 moves the readers).

They were one definition until the workflow stopped calling the module class:
the rule reads validation_psf members out of the per-exposure tars, keeps their
native dtypes and reconciles its output, none of which the class does or should
do. Two implementations is the right answer for the behaviour; two COLUMN LISTS
is not, and nothing else would notice them diverging — a column added to one
writer would simply be absent from the other's product, discovered by whoever
next tried to compute rho statistics from the wrong one.

Hence this module, which asserts the one thing they must share. It does NOT
assert the dtypes: the whole point of the hdf5 writer is that they differ (the
FITS one widens every float to 1D). Only the names, and their order.

Deliberately import-light on the workflow side: merge_star_cat.py pulls in h5py
and astropy, which the class does too, so a container-free run is not on offer
here and is not worth contorting for.
"""

import importlib.util
import sys
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = REPO_ROOT / "workflow" / "scripts"
SCRIPT = SCRIPTS / "merge_star_cat.py"


def _load_workflow_merge():
    """Import the rule's script by path — ``workflow/scripts`` is not a package.

    Its own imports (build_index, hdf5_reconcile, persist_exp) are siblings it
    reaches through ``sys.path[0]``, which is how the rule invokes it, so the
    directory goes on the path here too.
    """
    assert SCRIPT.exists(), f"{SCRIPT} not found; the rule calls it by path"
    sys.path.insert(0, str(SCRIPTS))
    try:
        spec = importlib.util.spec_from_file_location("_merge_star_cat", SCRIPT)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    finally:
        sys.path.remove(str(SCRIPTS))
    return module


@pytest.fixture(scope="module")
def writers():
    """The two column lists, each in the order its writer emits them."""
    h5py = pytest.importorskip("h5py")           # noqa: F841 - workflow dep
    pytest.importorskip("astropy")
    workflow = _load_workflow_merge()
    from shapepipe.modules.merge_starcat_package.merge_starcat import (
        MergeStarCatPSFEX,
    )
    # The class carries (output name, source column) pairs plus its optional
    # set and appends CCD_NB last; the script carries output names throughout.
    module_columns = (
        tuple(out for out, _ in MergeStarCatPSFEX._COLUMNS)
        + tuple(out for out, _ in MergeStarCatPSFEX._OPTIONAL)
        + ("CCD_NB",)
    )
    return module_columns, tuple(workflow.ALL_COLUMNS)


def test_column_names_and_order_agree(writers):
    """Same names, same order — the schema both products promise."""
    module_columns, workflow_columns = writers
    assert workflow_columns == module_columns


def test_sixteen_columns(writers):
    """The count is itself the documented contract (README, config.yaml)."""
    module_columns, workflow_columns = writers
    assert len(module_columns) == 16
    assert len(workflow_columns) == 16


def test_ccd_nb_is_last(writers):
    """CCD_NB is appended per input file rather than read from one, in both."""
    module_columns, workflow_columns = writers
    assert module_columns[-1] == "CCD_NB"
    assert workflow_columns[-1] == "CCD_NB"
