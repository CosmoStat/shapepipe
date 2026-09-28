"""The standalone ``scripts/python/merge_final_cat.py`` drops catalogue-param
columns absent from its input catalogues.

``TILE_UNIQUE_ID`` is absent from final catalogues made before make_cat wrote
it; ``filter_available_columns`` makes that column (and any other
requested-but-absent one) optional for this manual tool, logging one line per
drop. The workflow's campaign merge (``final_cat_merge``, through
``create_final_cat.read_data``) does the opposite and raises on a missing
column; tests/unit/test_final_cat_merge_invariants.py holds it to that.
"""

import importlib.util
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
MERGE_SCRIPT = REPO_ROOT / "scripts" / "python" / "merge_final_cat.py"


def _load(script):
    """Import a script by path — ``scripts/python`` is not a package."""
    assert script.exists(), f"{script} not found; the rule calls it by path"
    spec = importlib.util.spec_from_file_location(f"_{script.stem}", script)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_present_columns_are_kept():
    m = _load(MERGE_SCRIPT)
    param_list = ["XWIN_WORLD", "TILE_ID"]
    assert m.filter_available_columns(
        param_list, ["XWIN_WORLD", "TILE_ID", "FLAGS"]
    ) == param_list


def test_missing_column_is_dropped_and_logged(capsys):
    m = _load(MERGE_SCRIPT)
    kept = m.filter_available_columns(
        ["XWIN_WORLD", "TILE_UNIQUE_ID"], ["XWIN_WORLD", "FLAGS"]
    )
    assert kept == ["XWIN_WORLD"]
    out = capsys.readouterr().out
    assert "TILE_UNIQUE_ID" in out


def test_empty_param_list_means_copy_all_columns():
    m = _load(MERGE_SCRIPT)
    assert m.filter_available_columns([], ["XWIN_WORLD"]) == []
