"""A zero-detection tile traverses the real tile modules without fake rows."""

import importlib
import json
from pathlib import Path
import shutil
import sqlite3
import subprocess
import sys

from astropy.io import fits
import h5py
import numpy as np
import pytest
from sqlitedict import SqliteDict

from shapepipe.modules.sextractor_package import sextractor_script as sx
from shapepipe.modules.sextractor_package import match_catalogue as mc
from shapepipe.modules.merge_headers_package.merge_headers import merge_headers
from shapepipe.modules.vignetmaker_package.vignetmaker import VignetMaker
from shapepipe.modules.psfex_interp_package.psfex_interp import PSFExInterpolator
from shapepipe.modules.ngmix_runner import ngmix_runner
from shapepipe.modules.merge_sep_cats_package.merge_sep_cats import MergeSep
from shapepipe.modules.make_cat_runner import make_cat_runner
from shapepipe.pipeline.config import CustomParser
from tests.module.test_sextractor_post_process import (
    _write_exposure_headers, _write_sex_ldac,
)
from tests.module.test_match_catalogue import _external, _grid
from tests.module.test_ngmix_empty_tile import _RecordingLogger, _runner_config


@pytest.mark.parametrize("epoch_slots", [None, 12])
def test_zero_object_tile_modules_and_publication(tmp_path, monkeypatch, epoch_slots):
    scripts = Path(__file__).resolve().parents[2] / "workflow/scripts"
    monkeypatch.syspath_prepend(str(scripts))
    complete = importlib.import_module("completeness")
    tile = "177.257"
    number = "-177-257"
    run = tmp_path / tile
    log = _RecordingLogger()

    def output(stage, runner):
        path = run / "output" / stage / runner / "output"
        path.mkdir(parents=True, exist_ok=True)
        return path

    mdir = run / "manifests"
    mdir.mkdir(parents=True)

    def manifest(stage, subdir=None, suffix=""):
        record, ok = complete.build_manifest(stage, run, tile, subdir)
        assert ok, record["failures"]
        (mdir / f"{stage}{suffix}.json").write_text(json.dumps(record))
        return record

    # Fir 177.257's one exposure, but synthetic zero-row SExtractor output.
    exp = "2360341"
    npy = tmp_path / f"headers-{exp}.npy"
    _write_exposure_headers(npy)
    merge_headers([[str(npy)]], str(tmp_path), tile_number=number)
    headers = str(tmp_path / f"log_exp_headers{number}.sqlite")
    detected = output("run_sp_tile_Sx", "sextractor_runner")
    cat = detected / f"sexcat{number}.fits"
    _write_sex_ldac(cat, [exp], np.empty((0, 2)))
    with fits.open(cat, mode="update") as hdus:
        cols = hdus["LDAC_OBJECTS"].columns + fits.ColDefs([
            fits.Column(name=name, format="D", array=np.empty(0))
            for name in ["X_IMAGE", "Y_IMAGE", "X_IMAGE_DBL", "Y_IMAGE_DBL"]
        ] + [fits.Column(name="VIGNET", format="9E", dim="(3,3)",
                         array=np.empty((0, 3, 3), np.float32))])
        hdus[2] = fits.BinTableHDU.from_columns(cols, name="LDAC_OBJECTS")
    seg = tmp_path / "seg.fits"
    fits.writeto(seg, np.zeros((10, 10), np.int32))
    sx.add_seg_vignet(str(cat), str(seg), w_log=log)
    with fits.open(cat) as hdus:
        assert hdus[2].data["SEG_VIGNET"].shape == (0, 3, 3)
        assert hdus[2].data["SEG_VIGNET"].dtype.kind == "i"
    x, y = _grid(30)
    ext = _external(tmp_path / "external.cat", range(1, 31), x, y)
    counts = mc.match_catalogue(str(cat), ext, w_log=log)
    sx.make_post_process(str(cat), headers, ["XWIN_WORLD", "YWIN_WORLD"],
                         ["0", "100", "0", "100"])
    assert mc.read_match_record(cat) == counts
    (detected / "background.fits").touch()
    assert manifest("tile_detect")["catalogue_match"] == counts
    with fits.open(cat) as hdus:
        assert len(hdus["EPOCH_0"].data) == 0
        assert len(hdus["LDAC_OBJECTS"].data) == 0

    # Real PSF interpolation and all four multi-epoch vignette stores.
    psfout = output("run_sp_tile_PiViVi", "psfex_interp_runner")
    psf = PSFExInterpolator(None, str(cat), str(psfout), number, log,
                            pos_params=["XWIN_WORLD", "YWIN_WORLD"])
    psf.process_me([str(tmp_path)], "psf", headers)
    viout = output("run_sp_tile_PiViVi", "vignetmaker_runner")
    vi = VignetMaker(str(cat), "SPHE", ["XWIN_WORLD", "YWIN_WORLD"],
                     str(viout), number, log)
    vi.process_me([[str(tmp_path)]], ["image", "weight", "flag"], headers, 1)
    bgout = viout
    vi = VignetMaker(str(cat), "SPHE", ["XWIN_WORLD", "YWIN_WORLD"],
                     str(bgout), number, log)
    vi.process_me([[str(tmp_path)]], ["exp_background", "exp_background_rms"], headers, 1)
    stores = [viout / f"{name}_vignet{number}.sqlite"
              for name in ["image", "weight", "flag"]]
    psfpath = psfout / f"galaxy_psf{number}.sqlite"
    bgpath = bgout / f"exp_background_vignet{number}.sqlite"
    for path in [*stores, psfpath, bgpath]:
        with SqliteDict(str(path)) as db:
            assert len(db) == 0
    monkeypatch.setenv("SP_PSF", "psfex")
    manifest("tile_vignets")

    # Run both chunks, merge their five empty shear HDUs, assemble HDF5.
    for k in [1, 2]:
        subdir = f"run_sp_tile_ngmix_Ng{k}u"
        ngout = output(subdir, "ngmix_runner")
        log.messages.clear()
        config = _runner_config()
        config.set("NGMIX_RUNNER", "BLEND_HANDLING", "uberseg")
        config.set("NGMIX_RUNNER", "ID_OBJ_MIN", str(k))
        config.set("NGMIX_RUNNER", "ID_OBJ_MAX", str(k))
        assert ngmix_runner(
            [str(cat), str(stores[0]), str(bgpath), str(psfpath),
             str(stores[1]), str(stores[2]), headers],
            {"output": str(ngout)}, number, config, "NGMIX_RUNNER", log,
        ) == (None, None)
        logs = ngout.parent / "logs"
        logs.mkdir()
        (logs / f"process{number}.log").write_text("\n".join(log.messages))
        manifest("tile_ngmix", subdir, f"_{k}")
        with fits.open(ngout / f"ngmix{number}.fits") as hdus:
            assert len(hdus) == 6
            assert all(len(h.data) == 0 for h in hdus[1:])
    merged = output("run_sp_tile_Ms", "merge_sep_cats_runner")
    first = run / "output/run_sp_tile_ngmix_Ng1u/ngmix_runner/output" / f"ngmix{number}.fits"
    MergeSep([str(first)], number, ["ngmix"], [".fits"], str(merged),
             2, "error", log).process()
    manifest("tile_merge_cats")
    finalout = output("run_sp_tile_Mc", "make_cat_runner")
    config = CustomParser()
    options = {"SHAPE_MEASUREMENT_TYPE": "ngmix", "SAVE_PSF_DATA": "True"}
    if epoch_slots is not None:
        options["N_EPOCH_SLOTS"] = str(epoch_slots)
    config.read_dict({"MAKE_CAT_RUNNER": options})
    assert make_cat_runner(
        [str(cat), str(psfpath), str(merged / f"ngmix{number}.fits")],
        {"output": str(finalout)}, number, config, "MAKE_CAT_RUNNER", log,
    ) == (None, None)
    final = finalout / f"final_cat{number}.hdf5"
    with h5py.File(final) as h5:
        assert "TILE_ID" in h5 and "NGMIX_MCAL_FLAGS" in h5
        assert len(h5) > 50
        assert all(d.shape[0] == 0 for d in h5.values())
    monkeypatch.setenv("NGMIX_N_CHUNKS", "2")
    made = manifest("tile_make_cat")
    assert made["catalogue_match"] == counts
    assert made["epoch_cuts"]["totals"] == dict(
        considered=0, masked_fraction=0, central_veto=0, failed=0, objects_emptied=0)
    # Same persistent products as the publication rule (whose commands are
    # exercised separately by test_match_record_publication).
    published = tmp_path / "products/tiles/17" / tile
    published.mkdir(parents=True)
    shutil.copyfile(final, published / f"final_cat-{tile}.hdf5")
    shutil.copyfile(mdir / "tile_make_cat.json", published / "tile_make_cat.json")
    assert json.loads((published / "tile_make_cat.json").read_text())["catalogue_match"]["n_sextractor"] == 0

    # Campaign merge uses the shipped reader and reconciler, also used by
    # sp_latest; the tile remains present as a zero-length structured dataset.
    tiles = tmp_path / "tiles.txt"
    tiles.write_text(tile + "\n")
    index = tmp_path / "index.sqlite"
    with sqlite3.connect(index) as db:
        db.execute("CREATE TABLE tiles(tile_id TEXT PRIMARY KEY, ra_dir TEXT, n_exp INTEGER)")
        db.execute("INSERT INTO tiles VALUES (?, '177', 1)", (tile,))
    with h5py.File(final) as h5:
        columns = list(h5)
    params = tmp_path / "final_cat.param"
    params.write_text("\n".join(columns) + "\n")
    campaign_cat = tmp_path / "products/final_cat_empty.hdf5"
    argv = [sys.executable, str(scripts / "merge_final_cat.py"),
            "--products-dir", str(tmp_path / "products"), "--tile-list", str(tiles),
            "--index-db", str(index), "--output", str(campaign_cat),
            "--campaign", "empty", "--param-file", str(params)]
    merged = subprocess.run(argv, capture_output=True, text=True)
    assert merged.returncode == 0, merged.stdout + merged.stderr
    with h5py.File(campaign_cat) as h5:
        data = h5["patches/empty/177.257"]
        assert data.shape == (0,)
        assert set(data.dtype.names) == set(columns)
        assert h5.attrs["n_tiles"] == 1
    # An unchanged campaign must also pass the incremental-merge path.
    rerun = subprocess.run(argv, capture_output=True, text=True)
    assert rerun.returncode == 0, rerun.stderr
    assert "unchanged" in rerun.stdout
