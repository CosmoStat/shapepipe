"""Explicit numerical-frame contract for categorical registration."""

from pathlib import Path

import numpy as np
import numpy.testing as npt
import pytest
from astropy.coordinates import FK5, SkyCoord
from astropy.io import fits
from astropy.wcs import WCS

from shapepipe.modules.vignetmaker_package.vignetmaker import (
    normalise_registration_wcs,
    register_categories,
)

pytestmark = pytest.mark.decision("preparation.categorical_registration")


def _icrs():
    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    wcs.wcs.radesys = "ICRS"
    wcs.wcs.crval = [210.0, 63.5]
    wcs.wcs.crpix = [51.0, 51.0]
    wcs.wcs.cdelt = [-0.18577 / 3600, 0.18577 / 3600]
    return wcs


def _numeric_maps(seg, tile_wcs, epoch_wcs, centre, shape, number):
    """Numerical lookup with no frame reinterpretation or conversion."""
    rr, cc = np.indices(shape)
    xy = np.column_stack(
        ((cc + centre[1] - shape[1] // 2).ravel(),
         (rr + centre[0] - shape[0] // 2).ravel())
    )
    world = epoch_wcs.all_pix2world(xy, 0)
    tile = tile_wcs.all_world2pix(world, 0)
    x, y = np.floor(tile + .5).astype(int).T
    inside = (x >= 0) & (x < seg.shape[1]) & (y >= 0) & (y < seg.shape[0])
    expected = np.zeros(len(x), np.int32)
    expected[inside] = seg[y[inside], x[inside]]
    expected = expected.reshape(shape)
    return dict(SEG_VIGNET=expected, OFF_TILE=~inside.reshape(shape),
                NEIGHBOUR=(expected != 0) & (expected != number))


def test_real_gaia_fk5_header_preserves_numerical_lookup():
    path = Path(__file__).parents[1] / "data/registration_188307.npz"
    with np.load(path) as data:
        tile_header = fits.Header.fromstring(str(data["tile_header"]))
        epoch_header = fits.Header.fromstring(str(data["epoch_header"]))
        assert tile_header["RADESYS"] == "ICRS"
        assert epoch_header["RADECSYS"] == "FK5"
        assert epoch_header["EQUINOX"] == 2000.0
        assert epoch_header["ASTCAT"] == "test.gaia"
        tw, ew = WCS(tile_header), WCS(epoch_header)
        expected = _numeric_maps(
            data["segmentation"], tw, ew, data["int_pos"], (51, 51), 5422
        )
        got = register_categories(
            data["segmentation"], tw, ew, data["int_pos"], (51, 51), 5422,
            epoch_header=epoch_header,
        )
        for key in expected:
            npt.assert_array_equal(got[key], expected[key])
        npt.assert_array_equal(got["SEG_VIGNET"], data["expected_seg"])
        common = normalise_registration_wcs(tw, ew, epoch_header)
        assert common is not ew
        assert common.wcs.radesys == "ICRS"
        assert np.isnan(common.wcs.equinox)
        assert ew.wcs.radesys == "FK5" and ew.wcs.equinox == 2000.0
        assert epoch_header["RADECSYS"] == "FK5"
        pixels = np.array([[1., 2.], [1000.2, 1700.7]])
        npt.assert_array_equal(
            common.all_pix2world(pixels, 0), ew.all_pix2world(pixels, 0)
        )
        # A WCS alone cannot establish Gaia provenance.
        with pytest.raises(ValueError, match="equivalent celestial frames"):
            register_categories(
                data["segmentation"], tw, ew, data["int_pos"], (51, 51), 5422
            )


def test_genuinely_mixed_fk5_icrs_coordinates_raise():
    tw = _icrs()
    ew = tw.deepcopy()
    fk5 = SkyCoord(*tw.wcs.crval, unit="deg", frame="icrs").transform_to(
        FK5(equinox="J2000")
    )
    ew.wcs.crval = [fk5.ra.deg, fk5.dec.deg]
    ew.wcs.radesys = "FK5"
    ew.wcs.equinox = 2000.0
    seg = np.arange(101 * 101, dtype=np.int32).reshape(101, 101)
    with pytest.raises(ValueError, match="equivalent celestial frames"):
        register_categories(
            seg, tw, ew, [50, 50], (11, 11), 5100,
            epoch_header={"ASTCAT": "2MASS"},
        )


def test_icrs_control_matches_numerical_lookup():
    tw, ew = _icrs(), _icrs()
    ew.wcs.crpix += [.51, .51]
    seg = np.arange(101 * 101, dtype=np.int32).reshape(101, 101)
    expected = _numeric_maps(seg, tw, ew, [50, 50], (11, 11), 5100)
    assert normalise_registration_wcs(tw, ew) is ew
    got = register_categories(seg, tw, ew, [50, 50], (11, 11), 5100)
    for key in expected:
        npt.assert_array_equal(got[key], expected[key])


@pytest.mark.parametrize("astcat", ["test.gaia", "GAIA", "Gaia_DR3"])
def test_gaia_catalogue_token_normalises_fk5_j2000(astcat):
    tw, ew = _icrs(), _icrs()
    ew.wcs.radesys = "FK5"
    ew.wcs.equinox = 2000.0
    got = normalise_registration_wcs(tw, ew, {"ASTCAT": astcat})
    assert got.wcs.radesys == "ICRS" and np.isnan(got.wcs.equinox)


@pytest.mark.parametrize("astcat", ["", "2MASS", "notgaia", "gaia2"])
def test_non_gaia_catalogue_cannot_normalise_mismatched_frames(astcat):
    tw, ew = _icrs(), _icrs()
    ew.wcs.radesys = "FK5"
    ew.wcs.equinox = 2000.0
    with pytest.raises(ValueError, match="equivalent celestial frames"):
        normalise_registration_wcs(tw, ew, {"ASTCAT": astcat})


@pytest.mark.parametrize("frame,equinox", [("FK5", 2050.), ("FK4", 1950.)])
def test_other_mismatched_frames_raise_even_with_gaia(frame, equinox):
    tw, ew = _icrs(), _icrs()
    ew.wcs.radesys, ew.wcs.equinox = frame, equinox
    with pytest.raises(ValueError, match="equivalent celestial frames"):
        normalise_registration_wcs(tw, ew, {"ASTCAT": "test.gaia"})


def test_matching_fk5_frames_need_no_conversion():
    tw, ew = _icrs(), _icrs()
    for wcs in [tw, ew]:
        wcs.wcs.radesys = "FK5"
        wcs.wcs.equinox = 2000.0
    assert normalise_registration_wcs(tw, ew, {"ASTCAT": "2MASS"}) is ew
