"""``merge_final_cat.py`` and ``create_final_cat.py`` drop catalogue-param
columns absent from a given input catalogue.

``TILE_UNIQUE_ID`` is only present in per-tile catalogues produced with
``tile_detection: unions_catalogue``; a run using ``tile_detection:
sextractor`` never has it. Listing it in ``final_cat.param`` must not break
the merge for such runs: ``merge_final_cat.filter_available_columns`` and
``create_final_cat.read_data`` are what make that column (and any other
requested-but-absent column) optional, logging one line per drop instead of
failing partway through the merge.
"""

import importlib.util
from pathlib import Path

import numpy as np
from astropy.io import fits

REPO_ROOT = Path(__file__).resolve().parents[2]
MERGE_SCRIPT = REPO_ROOT / "scripts" / "python" / "merge_final_cat.py"
CREATE_SCRIPT = REPO_ROOT / "scripts" / "python" / "create_final_cat.py"


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


def _write_cat(path, columns):
    cols = [
        fits.Column(name=name, format="K", array=np.array(values, dtype=np.int64))
        for name, values in columns.items()
    ]
    hdu = fits.BinTableHDU.from_columns(cols)
    fits.HDUList([fits.PrimaryHDU(), hdu]).writeto(path, overwrite=True)


def test_create_final_cat_read_data_skips_missing_column(tmp_path, capsys):
    m = _load(CREATE_SCRIPT)

    cat_path = tmp_path / "final_cat-100-100.fits"
    _write_cat(cat_path, {"XWIN_WORLD": [1, 2], "TILE_ID": [100, 100]})

    params = {"hdu_num": 1, "param_list": ["XWIN_WORLD", "TILE_UNIQUE_ID"]}
    extracted_data, dtype, param_list = m.read_data(str(cat_path), params)

    assert param_list == ["XWIN_WORLD"]
    assert set(extracted_data) == {"XWIN_WORLD"}
    assert "TILE_UNIQUE_ID" not in dtype.names
    assert "TILE_UNIQUE_ID" in capsys.readouterr().out


def test_create_final_cat_read_data_keeps_present_column(tmp_path):
    m = _load(CREATE_SCRIPT)

    cat_path = tmp_path / "final_cat-100-100.fits"
    _write_cat(
        cat_path,
        {"XWIN_WORLD": [1, 2], "TILE_UNIQUE_ID": [100000001, 100000002]},
    )
    params = {"hdu_num": 1, "param_list": ["XWIN_WORLD", "TILE_UNIQUE_ID"]}
    extracted_data, dtype, param_list = m.read_data(str(cat_path), params)

    assert param_list == ["XWIN_WORLD", "TILE_UNIQUE_ID"]
    assert set(dtype.names) == {"XWIN_WORLD", "TILE_UNIQUE_ID"}
