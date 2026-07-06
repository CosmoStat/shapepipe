"""UNIT TESTS FOR VIGNETMAKER."""

from types import SimpleNamespace

from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
import numpy.testing as npt

from shapepipe.modules.vignetmaker_package.vignetmaker import (
    VignetMaker,
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

    stamps = maker._get_stamp(
        str(image_path),
        np.array([[10.0, 10.0], [14.2, 15.8]]),
        rad=2,
    )

    assert stamps.shape == (2, 5, 5)


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
