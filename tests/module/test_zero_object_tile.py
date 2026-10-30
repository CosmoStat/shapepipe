"""A zero-detection tile traverses the real tile modules without fake rows."""

import importlib
import json
from pathlib import Path
import shutil
import shlex
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
from shapepipe.modules.vignetmaker_runner import vignetmaker_runner
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
from tests.module.test_categorical_registration import _wcs


def test_zero_detection_check_image(tmp_path):
    """Real SExtractor writes a full-size all-zero segmentation even with no rows."""
    root = Path(__file__).resolve().parents[2] / "workflow/config/cfis"
    image = tmp_path / "image.fits"
    wcs = _wcs()
    fits.writeto(image, np.zeros((64, 64), np.float32), wcs.to_header())
    param = tmp_path / "empty.param"
    param.write_text("NUMBER\nX_IMAGE\nY_IMAGE\nVIGNET(3,3)\n")
    augmented = sx.seg_vignet_param_file(str(param), str(tmp_path / "seg.param"))
    caller = sx.SExtractorCaller(
        [str(image)], str(tmp_path), "-empty", str(root / "default_tile.sex"),
        augmented, str(root / "gauss_3.0_7x7.conv"),
        False, False, False, False, False, False, False,
        check_image=["BACKGROUND", "SEGMENTATION"],
    )
    result = subprocess.run(
        shlex.split(caller.make_command_line("source-extractor")),
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    with fits.open(caller.path_output_file) as hdus:
        assert len(hdus["LDAC_OBJECTS"].data) == 0
    segmentation = caller.check_paths["SEGMENTATION"]
    with fits.open(segmentation) as hdus:
        assert hdus[0].data.shape == (64, 64)
        assert hdus[0].data.dtype.kind == "i"
        assert not hdus[0].data.any()
        assert hdus[0].header["CTYPE1"] == "RA---TAN"
        assert hdus[0].header["CRPIX1"] == wcs.wcs.crpix[0]
    sx.add_seg_vignet(caller.path_output_file, segmentation)
    with fits.open(caller.path_output_file) as hdus:
        assert hdus["LDAC_OBJECTS"].data["SEG_LABEL"].shape == (0,)
        assert hdus["LDAC_OBJECTS"].data["SEG_VIGNET"].shape == (0, 3, 3)


@pytest.mark.parametrize("epoch_slots", [None, 12])
def test_zero_object_tile_modules_and_publication(
        tmp_path, monkeypatch, epoch_slots, *, tile="177.257", n_chunks=2,
        products_dir=None, exposure="2360341"):
    scripts = Path(__file__).resolve().parents[2] / "workflow/scripts"
    monkeypatch.syspath_prepend(str(scripts))
    complete = importlib.import_module("completeness")
    number = "-" + tile.replace(".", "-")
    products_dir = products_dir or tmp_path / "products"
    run = tmp_path / tile
    log = _RecordingLogger()
    log.debug = log.info  # Exposure-directory discovery also emits debug logs.

    def output(stage, runner):
        path = run / "output" / stage / runner / "output"
        path.mkdir(parents=True, exist_ok=True)
        return path

    mdir = run / "manifests"
    mdir.mkdir(parents=True, exist_ok=True)

    def manifest(stage, subdir=None, suffix=""):
        record, ok = complete.build_manifest(stage, run, tile, subdir)
        assert ok, record["failures"]
        (mdir / f"{stage}{suffix}.json").write_text(json.dumps(record))
        return record

    # Fir 177.257's one exposure, but synthetic zero-row SExtractor output.
    exp = exposure
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
    # Zero detections still produce a full-size, all-sky SEGMENTATION check
    # image, not an absent file or a zero-axis FITS image. SExtractor writes
    # it before the join/post-processing (see test_zero_detection_check_image).
    seg = detected / f"segmentation{number}.fits"
    fits.writeto(seg, np.zeros((10, 10), np.int32), _wcs().to_header())
    sx.add_seg_vignet(str(cat), str(seg), w_log=log)
    with fits.open(cat) as hdus:
        assert hdus[2].data["SEG_VIGNET"].shape == (0, 3, 3)
        assert hdus[2].data["SEG_VIGNET"].dtype.kind == "i"
        assert hdus[2].data["SEG_VIGNET"].dtype.itemsize == 4
        labels = hdus[2].data["SEG_LABEL"]
        assert labels.shape == (0,)
        assert labels.dtype.kind == "i" and labels.dtype.itemsize == 4
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
        assert hdus["LDAC_OBJECTS"].data["SEG_LABEL"].shape == (0,)
        assert hdus["LDAC_OBJECTS"].data["SEG_VIGNET"].shape == (0, 3, 3)

    # Real PSF interpolation and the four-input multi-epoch runners.
    psfout = output("run_sp_tile_PiViVi", "psfex_interp_runner")
    psf = PSFExInterpolator(None, str(cat), str(psfout), number, log,
                            pos_params=["XWIN_WORLD", "YWIN_WORLD"])
    psf.process_me([str(tmp_path)], "psf", headers)
    viout = output("run_sp_tile_PiViVi", "vignetmaker_runner")
    exp_root = tmp_path / "exp"
    exp_numbers = tmp_path / f"exp_numbers{number}.txt"
    exp_numbers.write_text(exp + "\n")
    # Directory discovery is real, but zero rows must not read CCD images.
    (exp_root / exp[:2] / exp / "output/run_sp_exp/split_exp_runner/output").mkdir(
        parents=True)
    config = CustomParser()
    config.read_dict({"VIGNETMAKER_RUNNER": dict(
        MASKING="False", STAMP_SIZE="3", COORD="SPHE", MODE="MULTI-EPOCH",
        POSITION_PARAMS="XWIN_WORLD,YWIN_WORLD", ME_IMAGE_EXP_DIR=str(exp_root),
        ME_IMAGE_EXP_RUNNERS="split_exp_runner",
    )})
    for patterns in ["image,weight,flag", "exp_background,exp_background_rms"]:
        config.set("VIGNETMAKER_RUNNER", "ME_IMAGE_PATTERN", patterns)
        assert vignetmaker_runner(
            [str(cat), headers, str(exp_numbers), str(seg)],
            {"output": str(viout)}, number, config, "VIGNETMAKER_RUNNER", log,
        ) == (None, None)
    bgout = viout
    stores = [viout / f"{name}_vignet{number}.sqlite"
              for name in ["image", "weight", "flag"]]
    psfpath = psfout / f"galaxy_psf{number}.sqlite"
    bgpath = bgout / f"exp_background_vignet{number}.sqlite"
    for path in [*stores, psfpath, bgpath,
                 bgout / f"exp_background_rms_vignet{number}.sqlite"]:
        assert path.is_file()
        with SqliteDict(str(path), flag="r") as db:
            # Stores are keyed by real object NUMBERs: no fake row or schema
            # sentinel. Category shapes/dtypes on actual image epochs are
            # checked by test_vignetmaker_writes_joined_categories_with_image_grid.
            assert len(db) == 0
    monkeypatch.setenv("SP_PSF", "psfex")
    manifest("tile_vignets")

    # Run both chunks, merge their five empty shear HDUs, assemble HDF5.
    for k in range(1, n_chunks + 1):
        subdir = f"run_sp_tile_ngmix_Ng{k}u"
        ngout = output(subdir, "ngmix_runner")
        log.messages.clear()
        config = _runner_config()
        config.set("NGMIX_RUNNER", "BLEND_HANDLING", "uberseg")
        # ngmix_range's bounded empty interval also works on registered stores.
        config.set("NGMIX_RUNNER", "ID_OBJ_MIN", "2")
        config.set("NGMIX_RUNNER", "ID_OBJ_MAX", "1")
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
             n_chunks, "error", log).process()
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
        assert not {"SEG_LABEL", "SEG_VIGNET", "VIGNET"}.intersection(h5)
        assert len(h5) > 50
        assert all(d.shape[0] == 0 for d in h5.values())
    monkeypatch.setenv("NGMIX_N_CHUNKS", str(n_chunks))
    made = manifest("tile_make_cat")
    assert made["catalogue_match"] == counts
    assert made["epoch_cuts"]["totals"] == dict(
        considered=0, masked_fraction=0, central_veto=0, failed=0, objects_emptied=0)
    # Same persistent products as the publication rule (whose commands are
    # exercised separately by test_match_record_publication).
    published = products_dir / "tiles" / tile[:2] / tile
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
        db.execute("INSERT INTO tiles VALUES (?, ?, 1)", (tile, tile.split(".")[0]))
    with h5py.File(final) as h5:
        columns = list(h5)
    params = tmp_path / "final_cat.param"
    params.write_text("\n".join(columns) + "\n")
    campaign_cat = products_dir / "final_cat_empty.hdf5"
    argv = [sys.executable, str(scripts / "merge_final_cat.py"),
            "--products-dir", str(products_dir), "--tile-list", str(tiles),
            "--index-db", str(index), "--output", str(campaign_cat),
            "--campaign", "empty", "--param-file", str(params)]
    merged = subprocess.run(argv, capture_output=True, text=True)
    assert merged.returncode == 0, merged.stdout + merged.stderr
    with h5py.File(campaign_cat) as h5:
        data = h5[f"patches/empty/{tile}"]
        assert data.shape == (0,)
        assert set(data.dtype.names) == set(columns)
        assert h5.attrs["n_tiles"] == 1
    # An unchanged campaign must also pass the incremental-merge path.
    rerun = subprocess.run(argv, capture_output=True, text=True)
    assert rerun.returncode == 0, rerun.stderr
    assert "unchanged" in rerun.stdout
