"""Empty edge tiles publish the same catalogue product as populated tiles."""

import subprocess
import sys

from astropy.io import fits
import numpy as np
import pytest
from sqlitedict import SqliteDict

from shapepipe.modules.ngmix_runner import ngmix_runner
from shapepipe.pipeline.config import CustomParser


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


class _RecordingLogger:
    """Minimal logger that records messages passed to warning/error/info."""

    def __init__(self):
        self.messages = []

    def info(self, msg, *_args, **_kwargs):
        self.messages.append(msg)

    def warning(self, msg, *_args, **_kwargs):
        self.messages.append(msg)

    def error(self, msg, *_args, **_kwargs):
        self.messages.append(msg)


def _runner_config():
    config = CustomParser()
    config.read_string(
        "[NGMIX_RUNNER]\n"
        "MAG_ZP = 30\n"
        "PIXEL_SCALE = .186\n"
        "ID_OBJ_MIN = -1\n"
        "ID_OBJ_MAX = -1\n"
    )
    return config


def test_psf_all_empty_guard_never_opens_galaxy_store(tmp_path):
    """The PSF-all-empty guard returns before the galaxy vignette store is
    opened at all -- the crash it guards against comes from unpickling that
    store's many-epoch arrays, so the check must fire without ever touching
    it, not just without acting on its contents.

    The galaxy ("image") store here is a corrupted, non-sqlite file:
    ``ngmix_runner`` must still succeed. Moving the PSF-empty guard after the
    image-vignet read -- or dropping it -- makes this test fail with a
    ``sqlite3`` error instead of the assertions below (checked by hand: with
    the two guards' order swapped, this test goes red).
    """
    ids = [1, 2, 3]
    psf_path = tmp_path / "galaxy_psf-001-001.sqlite"
    with SqliteDict(str(psf_path)) as db:
        for obj_id in ids:
            db[str(obj_id)] = {}
        db.commit()

    image_path = tmp_path / "image-001-001.sqlite"
    image_path.write_bytes(b"not a sqlite database")

    outputs = tmp_path / "output"
    outputs.mkdir()
    log = _RecordingLogger()

    input_file_list = [
        "tile_sexcat-001-001.fits",  # never opened: guard fires before Tile_cat
        str(image_path),             # galaxy store: corrupted, must stay unopened
        "exp_background-001-001.sqlite",
        str(psf_path),
        "weight-001-001.sqlite",
        "flag-001-001.sqlite",
        "log_exp_headers-001-001.sqlite",
    ]
    result = ngmix_runner(
        input_file_list,
        {"output": str(outputs)},
        "-001-001",
        _runner_config(),
        "NGMIX_RUNNER",
        log,
    )

    assert result == (None, None)
    with fits.open(outputs / "ngmix-001-001.fits") as hdul:
        assert {hdu.name for hdu in hdul[1:]} == {
            "NOSHEAR", "1P", "1M", "2P", "2M",
        }
        assert all(len(hdu.data) == 0 for hdu in hdul[1:])
    assert any(
        "all 3 objects failed the metacal fit (0 fitted)" in msg
        for msg in log.messages
    )
