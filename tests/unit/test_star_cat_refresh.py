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
        con.execute("CREATE TABLE tile_exposures(tile_id TEXT, exp_id TEXT)")
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
