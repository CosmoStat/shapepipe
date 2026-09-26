"""Which exposure-footprint records the workflow depends on, and how.

The footprint record (``exp_footprint.json``) is exp_footprint's declared
output, hanging off exp_persist's manifest and through it off the exposure's
scratch chain. Two rules name it — ``rule all`` and ``clean_exposure`` — through
``footprint_edge()``, which must drop it once the exposure's store is reclaimed:
a named record keeps the producer chain in the DAG, and any upstream params
change then reschedules the exposure's download, split and PSF fit.

The Snakefile helpers are lifted out by name and evaluated against a temporary
run root and products root, so these tests exercise what the helpers RETURN for
real files on disk.
"""

import re
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
SNAKEFILE = REPO_ROOT / "workflow" / "Snakefile"

EXP = "2243881"
HELPERS = ("exp_dir", "exp_manifest", "prod_exp_dir", "prod_exp_manifest",
           "tombstone", "exp_store_reclaimed", "footprint_edge")


def _snakefile_def(name):
    """The text of one top-level ``def`` in the Snakefile."""
    m = re.search(rf"^def {name}\(.*?(?=^\S)", SNAKEFILE.read_text(),
                  re.M | re.S)
    assert m, f"Snakefile no longer defines {name}()"
    return m.group(0)


@pytest.fixture
def roots(tmp_path):
    """The lifted helpers, bound to a temporary scratch and products root."""
    ns = {"Path": Path, "RUN_DIR": tmp_path / "run",
          "PRODUCTS_DIR": tmp_path / "products"}
    for name in HELPERS:
        exec(_snakefile_def(name), ns)
    return ns


def touch(path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("{}")


def test_live_store_depends_on_its_footprint(roots):
    """A live exposure is asked for its record: the thing to build, and what
    orders reclamation after the header read."""
    touch(roots["exp_manifest"](EXP, "exp_psf"))
    assert roots["footprint_edge"](EXP) == [
        roots["prod_exp_manifest"](EXP, "exp_footprint")]


def test_unbuilt_store_depends_on_its_footprint(roots):
    """A fresh campaign's exposure has no files yet and must still be asked."""
    assert roots["footprint_edge"](EXP) == [
        roots["prod_exp_manifest"](EXP, "exp_footprint")]


def test_cleaned_store_names_no_footprint(roots):
    """Tombstoned: the record is on the products root, and naming it would
    keep the exposure's producer chain in the DAG."""
    touch(roots["prod_exp_manifest"](EXP, "exp_persist"))
    touch(roots["prod_exp_manifest"](EXP, "exp_footprint"))
    touch(roots["tombstone"](EXP))
    assert roots["footprint_edge"](EXP) == []


def test_purged_store_names_no_footprint(roots):
    """Purged without a tombstone: persisted, and exp_psf's scratch manifest
    gone. The tombstone alone cannot see this case."""
    touch(roots["prod_exp_manifest"](EXP, "exp_persist"))
    touch(roots["prod_exp_manifest"](EXP, "exp_footprint"))
    assert roots["footprint_edge"](EXP) == []


def test_both_footprint_edges_use_footprint_edge():
    """``rule all`` and ``clean_exposure`` reach the record only through the
    helper, so the reclaimed-store cut cannot hold on one edge and not the
    other."""
    targets = _snakefile_def("footprint_targets")
    assert "footprint_edge(" in targets
    rules = (REPO_ROOT / "workflow" / "rules" / "exposure.smk").read_text()
    clean = re.search(r"^rule clean_exposure:.*?(?=^rule |\Z)", rules,
                      re.M | re.S).group(0)
    assert "footprint_edge(wc.exp)" in clean
    code = [line for src in (targets, clean) for line in src.splitlines()
            if not line.lstrip().startswith("#")]
    assert not any('"exp_footprint"' in line for line in code), code
