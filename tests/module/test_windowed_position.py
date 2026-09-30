"""UNIT TESTS FOR ``read_ext_sexcat_package.windowed_position``.

The windowed centroid read_ext_sexcat measures on UNIONS-catalogue tiles:
that it finds the centre of a Gaussian galaxy from a biased start, masks
neighbours as SExtractor does, falls back to the start where SExtractor
would, and reproduces SExtractor 2.25's own ``XWIN_IMAGE`` on real objects.
The fixture ``xwin_sextractor_186307.npz`` holds six objects of the r-band
tile CFIS.186.307 (DR5): cutouts of the image less SExtractor's background,
the SExtractor segmentation map relabelled to 1 for the object and 2, 3, ...
for neighbours, the starting ``X_IMAGE``, ``Y_IMAGE``, ``FLUX_RADIUS`` and
isophotal ellipse, and SExtractor's ``XWIN_IMAGE``, ``YWIN_IMAGE``,
``FLAGS_WIN``, all in cutout coordinates. SExtractor ran the committed
``default_tile.sex`` with the tile weight map.
"""

from pathlib import Path

import numpy as np
import pytest

from shapepipe.modules.read_ext_sexcat_package import windowed_position as wp

pytestmark = pytest.mark.decision("preparation.object_position_columns")

FIXTURE = Path(__file__).parent / "data" / "xwin_sextractor_186307.npz"


def _gaussian(shape, x_c, y_c, sigma, flux=1000.0):
    """Pixel-sampled circular Gaussian centred on 1-based ``(x_c, y_c)``."""
    y, x = np.mgrid[1:shape[0] + 1, 1:shape[1] + 1]
    r2 = (x - x_c) ** 2 + (y - y_c) ** 2
    return flux / (2 * np.pi * sigma**2) * np.exp(-r2 / (2 * sigma**2))


def test_gaussian_centre_found_from_a_truncated_barycentre():
    """A barycentre biased by a truncated isophote converges to the centre.

    The start is the barycentre of the pixels above a threshold on one side
    of the galaxy only (a truncated isophote), 0.8 sigma off; the window,
    sized by the true half-light radius, lands on the centre.
    """
    x_c, y_c, sigma = 30.37, 29.81, 2.0
    image = _gaussian((61, 61), x_c, y_c, sigma)
    y, x = np.mgrid[1:62, 1:62]
    iso = (image > 0.05 * image.max()) & (x >= x_c)
    x0, y0 = (np.sum(c[iso] * image[iso]) / image[iso].sum() for c in (x, y))
    assert np.hypot(x0 - x_c, y0 - y_c) > 1.0

    hlr = sigma * np.sqrt(2 * np.log(2))
    xw, yw, flags = wp.windowed_positions(image, [x0], [y0], [hlr])
    assert flags[0] == 0
    assert abs(xw[0] - x_c) < 1e-5 and abs(yw[0] - y_c) < 1e-5


def test_neighbour_footprint_is_masked_by_its_mirror():
    """A neighbour's footprint does not pull the centroid.

    A neighbour ten times brighter, 8 px away, captures the window when
    nothing masks it. With its footprint in the map, those pixels take the
    galaxy's mirror image and the centroid stays within 0.02 px of the
    isolated galaxy's (the rest is the neighbour's wings outside its
    footprint).
    """
    x_c, y_c, sigma = 30.0, 30.6, 2.0
    galaxy = _gaussian((61, 61), x_c, y_c, sigma)
    neighbour = _gaussian((61, 61), x_c + 8.0, y_c + 1.0, 1.5, flux=20000.0)
    y, x = np.mgrid[1:62, 1:62]
    seg = np.where(np.hypot(x - x_c - 8.0, y - y_c - 1.0) < 5.5, 2, 0)
    seg[np.hypot(x - x_c, y - y_c) < 3.0] = 1
    hlr = sigma * np.sqrt(2 * np.log(2))
    start = ([x_c + 0.5], [y_c], [hlr])

    alone = wp.windowed_positions(galaxy, *start)
    masked = wp.windowed_positions(galaxy + neighbour, *start, seg=seg,
                                   number=[1])
    unmasked = wp.windowed_positions(galaxy + neighbour, *start)
    assert np.hypot(masked[0][0] - alone[0][0],
                    masked[1][0] - alone[1][0]) < 0.02
    assert unmasked[0][0] - alone[0][0] > 1.0


def test_unmeasurable_and_negative_objects_keep_the_start():
    """Off-image, NaN and negative-flux objects keep the input position."""
    image = _gaussian((41, 41), 20.0, 20.0, 2.0)
    x0 = np.array([20.2, 60.0, np.nan])
    y0 = np.array([20.1, 20.0, 20.0])
    xw, yw, flags = wp.windowed_positions(image, x0, y0, [2.35, 2.35, 2.35])
    assert flags[0] == 0
    npt_equal = np.testing.assert_array_equal
    npt_equal(flags[1:], wp.FLAGS_WIN_UNMEASURED)
    npt_equal(xw[1:], x0[1:])

    xw, yw, flags = wp.windowed_positions(-image, [20.2], [20.1], [2.35])
    assert flags[0] & wp.FLAGS_WIN_NEGATIVE_FLUX
    assert (xw[0], yw[0]) == (20.2, 20.1)


def _fixture_objects():
    with np.load(FIXTURE) as data:
        names = sorted({k.rsplit("_", 1)[0] for k in data.files
                        if k.endswith("_image")})
        return {n: {k[len(n) + 1:]: data[k] for k in data.files
                    if k.startswith(n + "_")} for n in names}


@pytest.mark.parametrize("name", sorted(_fixture_objects()))
def test_matches_sextractor_on_real_objects(name):
    """SExtractor 2.25's XWIN_IMAGE and FLAGS_WIN, object by object.

    The objects are a bright and a faint isolated galaxy, two whose window
    reaches a neighbour's footprint (masked_neighbour*: unmasked, they are
    0.07 and 0.6 px off), one SExtractor falls back to the barycentre for
    (FLAGS_WIN 1) and one that has not converged after 16 steps. A fallback
    agrees only to the float32 rounding of the catalogue's X_IMAGE.
    """
    obj = _fixture_objects()[name]
    x0, y0, fr = obj["start"]
    xw, yw, flags = wp.windowed_positions(
        obj["image"], [x0], [y0], [fr], seg=obj["seg"], number=[1],
        cxx=obj["ellipse"][:1], cyy=obj["ellipse"][1:2],
        cxy=obj["ellipse"][2:],
    )
    expected = int(obj["flags_win"])
    assert flags[0] & 15 == expected
    tol = 1e-3 if expected else 1e-5
    assert abs(xw[0] - obj["expected"][0]) < tol
    assert abs(yw[0] - obj["expected"][1]) < tol
    if name == "not_converged":
        assert flags[0] & wp.FLAGS_WIN_NOT_CONVERGED


def test_background_recovers_a_smooth_sky():
    """A tilted sky plus sources comes back to well under the noise."""
    rng = np.random.default_rng(3)
    ny, nx = 1200, 1500
    y, x = np.mgrid[0:ny, 0:nx]
    sky = 5.0 + 1e-3 * x - 5e-4 * y
    image = sky + rng.normal(0, 1.0, (ny, nx))
    for _ in range(100):
        image += _gaussian((ny, nx), *rng.uniform(1, [nx, ny]), 2.0,
                           flux=rng.uniform(50, 500))
    back = wp.sextractor_background(image, back_size=256, filter_size=3)
    assert np.abs(back - sky)[100:-100, 100:-100].max() < 0.1
