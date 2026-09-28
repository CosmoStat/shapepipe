"""Static seam checks between HSM column producers and consumers.

Contracts ``psfex-validation-hsm-columns`` / ``psfex-starcat-columns-strict``
(psfex_interp ↔ merge_starcat) and ``psfex-me-shapes-columns`` /
``psf-epoch-slot-columns`` (psfex_interp ↔ make_cat). The producer side is
``_hsm_columns`` (``hsm-column-grammar``), which every psfex_interp writer
goes through; make_cat reads column names as string literals with no fallback,
collected from the function's AST, and merge_starcat reads its ``_COLUMNS``
table.
"""

import ast
import inspect
import textwrap

import numpy as np

from shapepipe.modules.make_cat_package.make_cat import SaveCatalogue
from shapepipe.modules.merge_starcat_package.merge_starcat import (
    MergeStarCatPSFEX,
)
from shapepipe.modules.psfex_interp_package.psfex_interp import (
    _HSM_ROW,
    _hsm_columns,
)


def _hsm_literals(func):
    """Set of ``HSM_*`` string literals in ``func``'s source."""
    tree = ast.parse(textwrap.dedent(inspect.getsource(func)))
    nodes = (n for n in ast.walk(tree) if isinstance(n, ast.Constant))
    return {
        n.value
        for n in nodes
        if isinstance(n.value, str) and n.value.startswith("HSM_")
    }


def _written(obj):
    return set(_hsm_columns(np.zeros((2, len(_HSM_ROW))), obj))


def test_merge_starcat_reads_exactly_what_psfex_validation_writes():
    """psfex-validation-hsm-columns == psfex-starcat-columns-strict.

    ``process`` reads every ``_COLUMNS`` source by name with no fallback (only
    ``_OPTIONAL`` is zero-filled), so the table is the read set.
    """
    written = _written("PSF") | _written("STAR")
    read = {
        col for _, col in MergeStarCatPSFEX._COLUMNS if col.startswith("HSM_")
    }
    assert not any(
        col.startswith("HSM_") for _, col in MergeStarCatPSFEX._OPTIONAL
    )
    assert written == read, {
        "written_not_read": sorted(written - read),
        "read_not_written": sorted(read - written),
    }


def test_make_cat_epoch_columns_come_from_psfex_me_shapes():
    """psf-epoch-slot-columns ⊆ psfex-me-shapes-columns."""
    produced = _written("PSF")
    consumed = _hsm_literals(SaveCatalogue._save_psf_data)
    assert consumed <= produced, sorted(consumed - produced)
