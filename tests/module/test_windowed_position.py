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


def test_negative_flux_flag_is_sextractors():
    """A centred negative source is FLAGS_WIN 4 exactly, as in SExtractor.

    SExtractor divides the (negative) second-moment sums by the negative
    flux, so the moments come out positive and only bit 4 is set.
    """
    image = -_gaussian((41, 41), 20.0, 20.0, 2.0)
    _, _, flags = wp.windowed_positions(image, [20.0], [20.0], [2.35])
    assert flags[0] == wp.FLAGS_WIN_NEGATIVE_FLUX


HLR = 2.0 * np.sqrt(2 * np.log(2))


def test_invalid_stripe_is_mirrored():
    """A no-data stripe through the window takes the mirror image.

    A Gaussian at (31, 31) with columns 33-35 holding no data: what the
    stripe holds is ignored, and a centroid at the centre stays there. Read
    as signal, the zero stripe drags it 0.6 px left. (From other starts
    SExtractor's integer mirror index can settle up to half a pixel off, at
    31.44 from (31.3, 30.8); this is SExtractor's behaviour too.)
    """
    image = _gaussian((61, 61), 31.0, 31.0, 2.0)
    valid = np.ones(image.shape, bool)
    valid[:, 32:35] = False
    zero = np.where(valid, image, 0.0)
    garbage = np.where(valid, image, 1e4)

    xw, yw, flags = wp.windowed_positions(zero, [31.0], [31.0], [HLR],
                                          valid=valid)
    assert flags[0] == 0
    assert abs(xw[0] - 31.0) < 1e-6 and abs(yw[0] - 31.0) < 1e-6
    for start in (([31.0], [31.0]), ([31.3], [30.8])):
        a = wp.windowed_positions(zero, *start, [HLR], valid=valid)
        b = wp.windowed_positions(garbage, *start, [HLR], valid=valid)
        np.testing.assert_array_equal(a[0], b[0])
        np.testing.assert_array_equal(a[1], b[1])
    xw, _, _ = wp.windowed_positions(zero, [31.0], [31.0], [HLR])
    assert xw[0] < 30.8


def test_invalid_and_off_image_mirrors_count_as_zero():
    """A mirror off the image, or itself invalid, contributes 0.

    Invalid pixels hold garbage; the result equals that of the same image
    with those pixels set to 0 and all valid, near the left edge (mirrors
    off the image) and with both a stripe and its mirror invalid.
    """
    rng = np.random.default_rng(5)
    for x_c, cols in ((4.0, np.r_[8:11]), (31.0, np.r_[26:29, 32:35])):
        image = _gaussian((61, 61), x_c, 31.0, 2.0)
        valid = np.ones(image.shape, bool)
        valid[:, cols] = False
        zeroed = np.where(valid, image, 0.0)
        garbage = np.where(valid, image, rng.uniform(1e3, 1e4, image.shape))
        start = ([x_c + 0.1], [31.2], [HLR])
        got = wp.windowed_positions(garbage, *start, valid=valid)
        want = wp.windowed_positions(zeroed, *start)
        np.testing.assert_allclose(got[:2], want[:2], rtol=0, atol=1e-9)
        assert got[2][0] == want[2][0]


def test_nan_outside_the_aperture_is_ignored():
    """A NaN that enters the stamp square after the first step changes nothing.

    From (30.2, 30.2) the square first spans 0-based columns 18-40; once
    the centre reaches (31, 31) it spans 19-41, and a NaN at column 41,
    outside the aperture, must not reach the sums.
    """
    image = _gaussian((61, 61), 31.0, 31.0, 2.0)
    clean = wp.windowed_positions(image, [30.2], [30.2], [HLR])
    image[30, 41] = np.nan
    got = wp.windowed_positions(image, [30.2], [30.2], [HLR])
    assert got[2][0] == 0
    np.testing.assert_allclose(got[:2], clean[:2], rtol=0, atol=1e-12)


def test_nan_inside_the_aperture_is_mirrored():
    """A NaN pixel on the galaxy is invalid and takes its mirror image."""
    image = _gaussian((61, 61), 31.0, 31.0, 2.0)
    image[31, 32] = np.nan
    xw, yw, flags = wp.windowed_positions(image, [30.6], [31.3], [HLR])
    assert flags[0] == 0
    assert abs(xw[0] - 31.0) < 1e-3 and abs(yw[0] - 31.0) < 1e-3


def test_wild_flux_radius_is_flagged_without_allocating():
    """FLUX_RADIUS <= 0 or NaN is an empty aperture; a huge one is skipped.

    DR6 carries FLUX_RADIUS down to -1.2e6; taken as a size, it built an
    84k-px box. Non-positive and NaN radii give FLAGS_WIN 4, radii whose
    window exceeds MAX_HALF give 16, both at the barycentre, and the
    measurable object alongside is unaffected.
    """
    import tracemalloc

    image = _gaussian((101, 101), 51.0, 51.0, 2.0)
    fr = [-12402.6, 0.0, np.nan, 1.18e6, HLR]
    x0 = [50.5, 50.5, 50.5, 50.5, 50.6]
    tracemalloc.start()
    xw, yw, flags = wp.windowed_positions(image, x0, [51.0] * 5, fr)
    peak = tracemalloc.get_traced_memory()[1]
    tracemalloc.stop()
    assert peak < 200e6
    np.testing.assert_array_equal(
        flags, [4, 4, 4, wp.FLAGS_WIN_UNMEASURED, 0])
    np.testing.assert_array_equal(xw[:4], x0[:4])
    assert abs(xw[4] - 51.0) < 1e-5


def test_box_is_clamped_to_a_small_image():
    """A window wider than the image measures the image, not a huge box."""
    image = _gaussian((9, 9), 5.0, 5.0, 1.0)
    xw, yw, flags = wp.windowed_positions(image, [5.2], [4.9], [60.0])
    assert flags[0] & wp.FLAGS_WIN_UNMEASURED == 0
    assert np.isfinite(xw[0])


def test_background_ignores_nan_pixels():
    """A NaN in a flat sky-3 mesh leaves the background at 3."""
    image = np.full((200, 200), 3.0)
    image[100, 100] = np.nan
    back = wp.sextractor_background(image, back_size=256, filter_size=3,
                                    good=image != 0)
    np.testing.assert_allclose(back, 3.0, atol=1e-6)


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
