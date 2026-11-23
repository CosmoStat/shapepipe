"""A PSF refit reaches the campaign star catalogue.

Real ``persist_exp.py`` and ``merge_star_cat.py``, run as the rules run them.
A refit that changes a value but no size rewrites the exposure's tar; its
manifest (the edge ``star_cat_merge`` waits on) must change with it, and the
next merge must refresh that exposure's dataset.
"""

import importlib.util
import json
import logging
import sqlite3
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np
from astropy.io import fits

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = REPO_ROOT / "workflow" / "scripts"
EXP = "2605805"


def _load(name):
    sys.path.insert(0, str(SCRIPTS))
    try:
        spec = importlib.util.spec_from_file_location(f"_{name}",
                                                      SCRIPTS / f"{name}.py")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    finally:
        sys.path.remove(str(SCRIPTS))
    return module


merge = _load("merge_star_cat")
persist = _load("persist_exp")


def _run(name, *args):
    run = subprocess.run([sys.executable, str(SCRIPTS / f"{name}.py"),
                          *map(str, args)], capture_output=True, text=True)
    assert run.returncode == 0, run.stderr
    return run.stdout


def _validation(path, value):
    """A PSFEx-shaped validation table: the rows live in HDU 2."""
    cols = [fits.Column(name=c, format="I" if "FLAG" in c else "E",
                        array=[value, value]) for c in merge.COLUMNS]
    path.parent.mkdir(parents=True, exist_ok=True)
    fits.HDUList([fits.PrimaryHDU(), fits.ImageHDU(),
                  fits.BinTableHDU.from_columns(cols)]).writeto(
        path, overwrite=True)


def test_validation_writer_emits_star_keys(monkeypatch, tmp_path):
    from shapepipe.modules.psfex_interp_package import psfex_interp as mod

    saved = {}

    class Output:
        def __init__(self, *args, **kwargs):
            pass

        def save_as_fits(self, data, **kwargs):
            saved.update(data)

    monkeypatch.setattr(mod.file_io, "FITSCatalogue", Output)
    monkeypatch.setattr(mod, "_hsm_columns", lambda *args: {})
    writer = mod.PSFExInterpolator.__new__(mod.PSFExInterpolator)
    writer._output_path_validation = str(tmp_path / "validation_psf")
    writer._img_number = "-02605805-12"
    writer._galcat_path = "input.fits"
    writer.psf_shapes = writer.star_shapes = None
    writer._write_output_validation({"NUMBER": np.array([9, 20])},
                                    {"SOURCE_NUMBER": [1, 2]})
    assert list(saved["EXP_ID"]) == ["02605805", "02605805"]
    assert list(saved["CCD_NB"]) == [12, 12]
    assert list(saved["SOURCE_ORDINAL"]) == [0, 1]


def test_star_key_unique(tmp_path):
    """Persist and merge two exposures, with two CCDs each, including leading zeros."""
    all_rows = []
    sources = []
    for exp in ("02605805", "2605806"):
        store = tmp_path / "scratch" / exp
        for ccd in (3, 12):
            src = (store / "output" / persist.RUN_NAME / "psfex_interp_runner"
                   / "output" / f"validation_psf-{exp}-{ccd}.fits")
            _validation(src, 0.25)
            sources.append([str(src)])
        dest = tmp_path / "exp" / exp[:2] / exp / "psf"
        manifest = dest.parent / "manifests" / "exp_persist.json"
        _run("persist_exp", "--exp-dir", store, "--exp", exp, "--dest", dest,
             "--manifest", manifest)
        all_rows.append(merge.read_exposure(exp, dest / f"{exp}.tar"))
    rows = np.concatenate(all_rows)
    assert {"EXP_ID", "CCD_NB", "SOURCE_ORDINAL"} <= set(rows.dtype.names)
    keys = [(exp.decode(), int(ccd), int(ordinal)) for exp, ccd, ordinal in
            zip(rows["EXP_ID"], rows["CCD_NB"], rows["SOURCE_ORDINAL"])]
    assert len(keys) == len(set(keys)) == 8
    assert set(keys) == {(exp, ccd, ordinal) for exp in ("02605805", "2605806")
                         for ccd in (3, 12) for ordinal in (0, 1)}

    from shapepipe.modules.merge_starcat_package.merge_starcat import MergeStarCatPSFEX
    MergeStarCatPSFEX(sources, str(tmp_path), logging.getLogger(__name__)).process()
    with fits.open(tmp_path / "full_starcat-0000000.fits") as hdul:
        flat = hdul[1].data
        assert len(flat) == len(keys)
        flat_keys = {(str(exp).strip(), int(ccd), int(ordinal)) for exp, ccd, ordinal
                     in zip(flat["EXP_ID"], flat["CCD_NB"], flat["SOURCE_ORDINAL"])}
    assert flat_keys == set(keys)


def test_refit_with_equal_sizes_refreshes_the_exposure(tmp_path):
    store = tmp_path / "scratch" / EXP
    src = (store / "output" / persist.RUN_NAME / "psfex_interp_runner"
           / "output" / f"validation_psf-{EXP}-3.fits")
    dest = tmp_path / "exp" / EXP[:2] / EXP / "psf"
    manifest = dest.parent / "manifests" / "exp_persist.json"
    persist_args = ("--exp-dir", store, "--exp", EXP, "--dest", dest,
                    "--manifest", manifest)

    tiles = tmp_path / "tiles.txt"
    tiles.write_text("210.282\n")
    db = tmp_path / "index.sqlite"
    with sqlite3.connect(db) as con:
        con.execute("CREATE TABLE tiles(tile_id TEXT PRIMARY KEY, "
                    "ra_dir TEXT, n_exp INTEGER)")
        con.execute("CREATE TABLE tile_exposures(tile_id TEXT, exp_id TEXT)")
        con.execute("INSERT INTO tiles VALUES ('210.282', '210', 1)")
        con.execute("INSERT INTO tile_exposures VALUES ('210.282', ?)", (EXP,))
    out = tmp_path / "stars.h5"
    merge_args = ("--products-dir", tmp_path, "--tile-list", tiles,
                  "--index-db", db, "--output", out, "--campaign", "t")

    _validation(src, 0.25)
    _run("persist_exp", *persist_args)
    _run("merge_star_cat", *merge_args)
    size, before = src.stat().st_size, manifest.read_bytes()

    _validation(src, 0.75)
    assert src.stat().st_size == size, "fixture must keep the size"
    _run("persist_exp", *persist_args)
    assert manifest.read_bytes() != before, (
        "the tar changed but the manifest star_cat_merge waits on did not")

    assert "1 refreshed" in _run("merge_star_cat", *merge_args)
    with h5py.File(out) as f:
        rows = f[f"exposures/{EXP}"][:]
    np.testing.assert_allclose(rows["X"], [0.75, 0.75])
    assert list(rows["EXP_ID"]) == [EXP.encode()] * 2
    assert list(rows["CCD_NB"]) == [3, 3]
    assert list(rows["SOURCE_ORDINAL"]) == [0, 1]
