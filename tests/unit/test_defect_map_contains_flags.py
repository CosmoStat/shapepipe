"""An exposure's defect fragment contains its flag image.

``rasterize_ccd`` (``workflow/scripts/defect_map_exp.py``) turns one CCD's flag
image and WCS into the healpix pixels its flags touch. The map it feeds is
CONSERVATIVE by contract (``defect-fragment-contains-flags``): a flagged CCD
pixel that lands in an unmasked healpix pixel is a hole in the footprint that
nothing downstream can see. So the invariant here is containment, checked
against the sky geometry itself and not against the function's own sampling:

  * every flagged pixel's centre AND its four corners — the vertices that bound
    its footprint — lie in masked healpix pixels, with the sky positions taken
    through astropy's 0-based ``pixel_to_world_values``, independent of the
    1-based ``all_pix2world`` path the rasterizer uses;
  * a clean flag image yields no pixels at all.

The fixture carries the three shapes a MegaCam flag image has: a ONE-PIXEL bad
column (the geometry that centre sampling erases), a saturated blob, and
isolated pixels, on the first and last row and column and on a lattice of hot
pixels spaced wider than a healpix pixel. The lattice is what makes the test
sharp: a lone pixel straddling a healpix boundary has no flagged neighbour
whose centre already masks the far side, so centre-only sampling misses it. The WCS is a rotated TAN
at 0.187 arcsec/pixel, MegaCam's scale, so a healpix pixel at the ladder's nside
covers ~74 CCD pixels, as on the sky. ``CHUNK`` is shrunk so the batching runs
several batches. The oversample is read from ``workflow/config.yaml``, the
value a campaign runs with.

Needs healpy, healsparse and astropy, so it runs inside the container and
skips outside.
"""

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest


pytestmark = pytest.mark.unions

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = REPO_ROOT / "workflow" / "scripts"
SCRIPT = SCRIPTS / "defect_map_exp.py"
CONFIG = REPO_ROOT / "workflow" / "config.yaml"

NY, NX = 320, 240                # a few hundred pixels, not a 2048 x 4612 chip
SCALE_DEG = 0.187 / 3600.0       # MegaCam
ROTATION_DEG = 23.0              # off-axis, so healpix boundaries cut obliquely
CRVAL = (150.3, 31.7)


def _load():
    """Import the rule's script by path — scripts/ is not a package."""
    assert SCRIPT.exists(), f"{SCRIPT} not found; the rule calls it by path"
    spec = importlib.util.spec_from_file_location("_defect_map_exp", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def raster():
    pytest.importorskip("healpy")
    pytest.importorskip("healsparse")
    pytest.importorskip("astropy")
    return _load()


@pytest.fixture(scope="module")
def ladder():
    """``(nside, oversample)`` as the campaign config sets them."""
    yaml = pytest.importorskip("yaml")
    block = yaml.safe_load(CONFIG.read_text())["defect_map"]
    return int(block["nside"]), int(block["oversample"])


def _wcs():
    from astropy.wcs import WCS

    theta = np.deg2rad(ROTATION_DEG)
    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    wcs.wcs.crval = list(CRVAL)
    wcs.wcs.crpix = [NX / 2 + 0.5, NY / 2 + 0.5]
    wcs.wcs.cd = SCALE_DEG * np.array([[-np.cos(theta), np.sin(theta)],
                                       [np.sin(theta), np.cos(theta)]])
    return wcs


def _flags():
    flags = np.zeros((NY, NX), dtype=np.int16)
    flags[:, 101] = 1                                   # one-pixel bad column
    yy, xx = np.mgrid[:NY, :NX]
    flags[(yy - 210) ** 2 + (xx - 60) ** 2 <= 7 ** 2] = 2   # saturated blob
    for row, col in [(0, 0), (0, NX - 1), (NY - 1, 0), (NY - 1, NX - 1),
                     (NY - 1, 150), (40, NX - 1), (77, 33)]:
        flags[row, col] = 8                             # isolated pixels
    # Hot pixels spaced wider than a healpix pixel (11 x 0.187" > 1.61"), so
    # one straddling a healpix boundary has no flagged neighbour whose centre
    # already masks the far side. This is what makes centre sampling fail.
    flags[5::11, 3::11] = 8
    return flags


def _write(tmp_path: Path, flags):
    """The two files ``rasterize_ccd`` reads: the flag split and the image
    split whose header alone carries the WCS."""
    from astropy.io import fits

    flag_path = tmp_path / "flag-2079612-0.fits"
    image_path = tmp_path / "image-2079612-0.fits"
    fits.PrimaryHDU(data=flags).writeto(flag_path)
    fits.PrimaryHDU(header=_wcs().to_header()).writeto(image_path)
    return flag_path, image_path


def _rasterize(raster, ladder, tmp_path, flags, monkeypatch):
    nside, oversample = ladder
    monkeypatch.setattr(raster, "CHUNK", 97)
    off_x, off_y = raster.offsets(oversample)
    flag_path, image_path = _write(tmp_path, flags)
    return raster.rasterize_ccd(flag_path, image_path, nside, off_x, off_y)


def _vertex_pixels(nside, flags):
    """``(healpix id per vertex, (row, col) per vertex)`` for the centre and the
    four corners of every flagged pixel, through astropy's 0-based convention."""
    import healpy as hp

    rows, cols = np.nonzero(flags)
    offsets = [(0.0, 0.0), (-0.5, -0.5), (-0.5, 0.5), (0.5, -0.5), (0.5, 0.5)]
    dx = np.array([o[0] for o in offsets])
    dy = np.array([o[1] for o in offsets])
    x = (cols[:, None] + dx[None, :]).ravel()
    y = (rows[:, None] + dy[None, :]).ravel()
    ra, dec = _wcs().pixel_to_world_values(x, y)
    ids = hp.ang2pix(nside, ra, dec, lonlat=True, nest=True)
    where = np.repeat(np.stack([rows, cols], axis=1), len(offsets), axis=0)
    return ids, where


def test_fragment_contains_every_flagged_pixel(raster, ladder, tmp_path,
                                               monkeypatch):
    """@sc defect-fragment-contains-flags: no flagged vertex lands unmasked."""
    nside, _ = ladder
    flags = _flags()
    masked = _rasterize(raster, ladder, tmp_path, flags, monkeypatch)

    ids, where = _vertex_pixels(nside, flags)
    missed = ~np.isin(ids, masked)
    assert not missed.any(), (
        f"{missed.sum()} of {missed.size} flagged-pixel vertices fall in "
        f"unmasked healpix pixels, at (row, col) "
        f"{sorted(set(map(tuple, where[missed].tolist())))[:10]}")


def test_clean_flag_image_masks_nothing(raster, ladder, tmp_path,
                                        monkeypatch):
    masked = _rasterize(raster, ladder, tmp_path,
                        np.zeros((NY, NX), dtype=np.int16), monkeypatch)
    assert masked.size == 0
