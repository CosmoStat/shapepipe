"""A PSF refit reaches the campaign star catalogue.

Real ``persist_exp.py`` and ``merge_star_cat.py``, run as the rules run them.
A refit that changes a value but no size rewrites the exposure's tar; its
manifest (the edge ``star_cat_merge`` waits on) must change with it, and the
next merge must refresh that exposure's dataset.
"""

import importlib.util
import json
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


def _mccd_validation(path):
    """An MCCD-shaped validation table: one exposure, HDU 1, moment vectors."""
    mom_psf = np.array([[0.01, -0.02, 1.5, 0], [0.03, 0.04, 2.0, 1]])
    mom_star = np.array([[0.02, -0.01, 1.6, 0], [0.05, 0.06, 2.1, 0]])
    cols = [
        fits.Column(name="PSF_MOM_LIST", format="4D", array=mom_psf),
        fits.Column(name="STAR_MOM_LIST", format="4D", array=mom_star),
        fits.Column(name="GLOB_POSITION_IMG_LIST", format="2D",
                    array=[[100.0, 200.0], [300.0, 400.0]]),
        fits.Column(name="CCD_ID_LIST", format="1D", array=[3.0, 17.0]),
        fits.Column(name="RA_LIST", format="1D", array=[10.0, 10.1]),
        fits.Column(name="DEC_LIST", format="1D", array=[30.0, 30.1]),
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    fits.HDUList([fits.PrimaryHDU(),
                  fits.BinTableHDU.from_columns(cols)]).writeto(
        path, overwrite=True)


def test_mccd_validation_merges_into_the_same_columns(tmp_path):
    """psf_model=mccd: one exposure-wide table becomes the PSFEx columns."""
    store = tmp_path / "scratch" / EXP
    src = (store / "output" / persist.RUN_NAME / "mccd_fit_val_runner"
           / "output" / f"validation_psf-{EXP}.fits")
    model = src.with_name(f"fitted_model-{EXP}.npy")
    dest = tmp_path / "exp" / EXP[:2] / EXP / "psf"
    manifest = dest.parent / "manifests" / "exp_persist.json"

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

    _mccd_validation(src)
    np.save(model, np.zeros(3))
    _run("persist_exp", "--exp-dir", store, "--exp", EXP, "--dest", dest,
         "--manifest", manifest, "--pattern", "mccd_model")
    packed = {f["name"]: f["product"]
              for f in json.loads(manifest.read_text())["files"]}
    assert packed == {f"validation_psf-{EXP}.fits": "psf_validation",
                      f"fitted_model-{EXP}.npy": "mccd_model"}

    _run("merge_star_cat", "--products-dir", tmp_path, "--tile-list", tiles,
         "--index-db", db, "--output", out, "--campaign", "t")
    with h5py.File(out) as f:
        rows = f[f"exposures/{EXP}"][:]
    assert rows.dtype.names == merge.ALL_COLUMNS
    np.testing.assert_allclose(rows["X"], [100.0, 300.0])
    np.testing.assert_allclose(rows["DEC"], [30.0, 30.1])
    np.testing.assert_allclose(rows["HSM_G2_PSF"], [-0.02, 0.04])
    np.testing.assert_allclose(rows["HSM_T_PSF"], 2 * np.array([1.5, 2.0])**2)
    np.testing.assert_allclose(rows["HSM_T_STAR"], 2 * np.array([1.6, 2.1])**2)
    np.testing.assert_array_equal(rows["HSM_FLAG_PSF"], [0, 1])
    np.testing.assert_array_equal(rows["CCD_NB"], [3, 17])
    assert np.isnan(rows["HSM_RHO4_STAR"]).all()
    np.testing.assert_array_equal(rows["ACCEPTED"], [0, 0])
