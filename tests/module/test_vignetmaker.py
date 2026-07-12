"""UNIT TESTS FOR VIGNETMAKER."""

from types import SimpleNamespace

from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
import numpy.testing as npt
import pytest

from shapepipe.modules.vignetmaker_package.vignetmaker import (
    VignetMaker,
    get_stamps,
    make_mask,
)


def _make_tan_header():

    header = fits.Header()
    header["NAXIS"] = 2
    header["NAXIS1"] = 21
    header["NAXIS2"] = 21
    header["CTYPE1"] = "RA---TAN"
    header["CTYPE2"] = "DEC--TAN"
    header["CRPIX1"] = 11.0
    header["CRPIX2"] = 11.0
    header["CRVAL1"] = 180.0
    header["CRVAL2"] = 30.0
    header["CD1_1"] = -1.0e-3
    header["CD1_2"] = 0.0
    header["CD2_1"] = 0.0
    header["CD2_2"] = 1.0e-3

    return header


def _write_sex_catalogue(path, **columns):

    def make_column(name, values):

        values = np.asarray(values)
        if values.ndim == 1:
            return fits.Column(name=name, format="D", array=values)
        return fits.Column(
            name=name,
            format=f"{np.prod(values.shape[1:])}D",
            dim=str(values.shape[1:]),
            array=values,
        )

    fits.HDUList(
        [
            fits.PrimaryHDU(),
            fits.BinTableHDU(),
            fits.BinTableHDU.from_columns(
                [make_column(name, values) for name, values in columns.items()]
            ),
        ]
    ).writeto(path)


def test_convert_pos_roundtrips_world_coordinates(tmp_path):

    image_path = tmp_path / "image.fits"
    header = _make_tan_header()
    fits.PrimaryHDU(data=np.zeros((21, 21)), header=header).writeto(image_path)

    pixels_xy = np.array([[8.0, 9.0], [11.0, 11.0], [14.0, 13.0]])
    worlds = WCS(header).all_pix2world(pixels_xy[:, [1, 0]], 1)[:, [1, 0]]
    maker = VignetMaker.__new__(VignetMaker)
    maker._pos = worlds

    converted = maker.convert_pos(str(image_path))

    npt.assert_allclose(converted, pixels_xy, atol=1e-6)


def test_get_stamp_preserves_count_and_stamp_shape(tmp_path):

    image_path = tmp_path / "image.fits"
    image = np.arange(25 * 25, dtype=np.float32).reshape(25, 25)
    fits.PrimaryHDU(data=image).writeto(image_path)

    maker = VignetMaker.__new__(VignetMaker)
    maker._w_log = SimpleNamespace(info=lambda message: None)

    positions = np.array([[10.0, 10.0], [14.2, 15.8]])
    stamps, int_pos, offsets = maker._get_stamp(
        str(image_path),
        positions,
        rad=2,
    )

    assert stamps.shape == (2, 5, 5)
    npt.assert_array_equal(int_pos, np.round(positions).astype(int))
    npt.assert_allclose(offsets, positions - np.round(positions))


def test_get_stamps_extracts_around_rounded_pixel():
    """In-bounds stamps equal a manual zero-padded slice; bookkeeping holds."""
    rng = np.random.default_rng(0)
    image = rng.normal(size=(40, 40)).astype(np.float32)
    rad = 3
    positions = np.array([[10.0, 12.0], [20.4, 5.6], [0.0, 39.0], [39.0, 0.0]])

    stamps, int_pos, offsets = get_stamps(image, positions, rad)

    assert stamps.shape == (4, 2 * rad + 1, 2 * rad + 1)
    npt.assert_array_equal(int_pos, np.round(positions).astype(int))
    npt.assert_allclose(offsets, positions - np.round(positions))

    padded = np.pad(image, ((rad, rad), (rad, rad)), mode="constant")
    for stamp, (row, col) in zip(stamps, int_pos):
        npt.assert_array_equal(
            stamp, padded[row : row + 2 * rad + 1, col : col + 2 * rad + 1]
        )


def test_get_stamps_zero_pads_edges():
    """A corner object is zero-padded on the out-of-frame side."""
    image = np.ones((10, 10))
    rad = 2

    stamps, _, _ = get_stamps(image, np.array([[0.0, 0.0]]), rad)

    expected = np.ones((5, 5))
    expected[:2, :] = 0.0  # rows above the top edge
    expected[:, :2] = 0.0  # cols left of the left edge
    npt.assert_array_equal(stamps[0], expected)


def test_get_stamps_offset_is_subpixel():
    """Offset is the catalog float minus its rounded pixel, in [-0.5, 0.5]."""
    image = np.zeros((20, 20))
    positions = np.array([[10.3, 10.0], [10.7, 9.4], [10.5, 10.5]])

    _, int_pos, offsets = get_stamps(image, positions, 2)

    npt.assert_allclose(offsets, positions - int_pos)
    assert np.all(np.abs(offsets) <= 0.5 + 1e-12)


@pytest.mark.parametrize(
    "positions",
    [
        np.array([[-1.0, 5.0]]),  # above the top edge
        np.array([[5.0, 10.0]]),  # past the right edge (shape is 10)
        np.array([[5.0, 5.0], [10.0, 5.0]]),  # one good, one out
    ],
)
def test_get_stamps_raises_out_of_bounds(positions):
    """Unlike sf_tools' modulo wrap, an out-of-bounds centre raises."""
    image = np.zeros((10, 10))
    with pytest.raises(ValueError, match="outside the image"):
        get_stamps(image, positions, 2)


@pytest.mark.parametrize("rad", [2, 3, 5])
def test_get_stamps_matches_fetchstamps(rad):
    """Bit-identical to sf_tools FetchStamps for in-bounds objects.

    The extraction is pure slicing of a zero-padded array, so it is
    value-independent: random data with positions spanning the interior and
    the padded edges exercises exactly what a real exposure would. Skipped
    once sf_tools is no longer installed (it is being dropped as a dep).
    """
    stamp_mod = pytest.importorskip("sf_tools.image.stamp")
    rng = np.random.default_rng(42)
    image = rng.normal(size=(60, 60))
    # In-bounds positions, including edges (padding) — but never out of
    # bounds, where FetchStamps wraps and get_stamps deliberately raises.
    positions = rng.uniform(0.0, 59.0, size=(50, 2))

    stamps, _, _ = get_stamps(image, positions, rad)

    fetch = stamp_mod.FetchStamps(image, int(rad))
    fetch.get_pixels(np.round(positions).astype(int))
    npt.assert_array_equal(stamps, fetch.scan())


def test_make_mask_replaces_only_sentinel_pixels(tmp_path):

    cat_path = tmp_path / "galaxies.fits"
    vignets = np.array(
        [
            [[1.0, -1.0e30], [3.0, -1.0e29]],
            [[-2.0e30, 5.0], [6.0, 7.0]],
        ]
    )
    _write_sex_catalogue(cat_path, VIGNET=vignets)

    masked = make_mask(str(cat_path), mask_value=-99.0)

    npt.assert_array_equal(
        masked,
        np.array(
            [
                [[1.0, -99.0], [3.0, -1.0e29]],
                [[-99.0, 5.0], [6.0, 7.0]],
            ]
        ),
    )
