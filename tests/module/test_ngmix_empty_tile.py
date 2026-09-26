"""Empty edge tiles publish the same catalogue product as populated tiles."""

import subprocess
import sys

from astropy.io import fits
import numpy as np
import pytest
from sqlitedict import SqliteDict


@pytest.mark.parametrize("empty_store", ["galaxy", "psf"])
def test_empty_tile_cli_writes_catalogue_and_health_log(tmp_path, empty_store):
    """Contract empty-tile-product through the shipped launcher.

    A successful exit without a catalogue fails the campaign's per-tile
    completeness check. Both empty-galaxy and empty-PSF stores must publish
    one catalogue with five empty metacal HDUs and log the zero-fitted run.
    """
    inputs = tmp_path / "input"
    outputs = tmp_path / "output"
    inputs.mkdir()
    outputs.mkdir()
    ids = np.arange(1, 4)
    objects = fits.BinTableHDU.from_columns(
        [
            fits.Column(name="NUMBER", format="J", array=ids),
            fits.Column(name="XWIN_WORLD", format="D", array=np.zeros(3)),
            fits.Column(name="YWIN_WORLD", format="D", array=np.zeros(3)),
        ],
        name="LDAC_OBJECTS",
    )
    imhead = fits.BinTableHDU.from_columns(
        [fits.Column(name="Field Header Card", format="1A", array=["x"])],
        name="LDAC_IMHEAD",
    )
    fits.HDUList([fits.PrimaryHDU(), imhead, objects]).writeto(
        inputs / "tile_sexcat-001-001.fits"
    )
    epoch = {"exp-1": {"VIGNET": np.ones((5, 5)), "OFFSET": [0.0, 0.0]}}
    stores = {
        "image": "empty" if empty_store == "galaxy" else epoch,
        "galaxy_psf": {} if empty_store == "psf" else epoch,
        "exp_background": {},
        "weight": {},
        "flag": {},
        "log_exp_headers": {},
    }
    for name, value in stores.items():
        with SqliteDict(str(inputs / f"{name}-001-001.sqlite")) as db:
            for obj_id in ids:
                db[str(obj_id)] = value
            db.commit()

    config = tmp_path / "config.ini"
    config.write_text(f"""[DEFAULT]
RUN_NAME = run_empty
RUN_DATETIME = False
VERBOSE = False
[EXECUTION]
MODULE = ngmix_runner
MODE = smp
[FILE]
INPUT_DIR = {inputs}
OUTPUT_DIR = {outputs}
NUMBERING_SCHEME = -000-000
[JOB]
SMP_BATCH_SIZE = 1
TIMEOUT = 00:02:00
[NGMIX_RUNNER]
MAG_ZP = 30
PIXEL_SCALE = .186
ID_OBJ_MIN = -1
ID_OBJ_MAX = -1
""")
    run = subprocess.run(
        [sys.executable, "-m", "shapepipe.shapepipe_run", "-c", str(config)],
        cwd=tmp_path, capture_output=True, text=True, timeout=120,
    )
    assert run.returncode == 0, run.stdout + run.stderr
    product_dir = outputs / "run_empty" / "ngmix_runner" / "output"
    products = list(product_dir.iterdir())
    assert [path.name for path in products] == ["ngmix-001-001.fits"]
    with fits.open(products[0]) as hdul:
        assert {hdu.name for hdu in hdul[1:]} == {
            "NOSHEAR", "1P", "1M", "2P", "2M",
        }
        assert all(len(hdu.data) == 0 for hdu in hdul[1:])
    log = "\n".join(path.read_text() for path in outputs.rglob("process*.log"))
    assert log.count("all 3 objects failed the metacal fit (0 fitted)") == 1
