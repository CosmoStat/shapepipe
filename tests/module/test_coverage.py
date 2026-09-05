"""UNIT / PROPERTY TESTS FOR THE COVERAGE-MASK GEOMETRY.

Covers the pure logic the coverage nexp mask rests on, on both sides of the
``exp_footprint`` -> ``coverage_map`` handover: per-CCD image-shape resolution
from plain and fpack-compressed headers, the sky-corner projection, and then
the accumulated exposure-count (nexp) contract, the RA-wrap and pole guards,
the power-of-two ``nside`` validation and the median smoothing of a finished
map. Everything here takes arrays or headers in memory; the way a record
reaches these functions is ``tests/unit/test_exp_footprint.py``'s subject, and
the heavy paths (production-resolution map building, plotting) are exercised
end to end by the pipeline.
"""

import healsparse as hsp
import hpgeom as hpg
import numpy as np
import numpy.testing as npt
import pytest
from astropy import wcs
from astropy.io.fits import Header

from shapepipe.utilities.ccd_footprint import _ccd_corners, _image_shape
from shapepipe.utilities.coverage_map_builder import (
    build_map,
    check_nside,
    median_filter,
    unwrap_ra,
)


def _tan_wcs(crval_ra, crval_dec=0.0, nx=2080, ny=4612):
    """Build a minimal 2-D TAN WCS with a populated pixel shape."""
    w = wcs.WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [crval_ra, crval_dec]
    w.wcs.crpix = [nx / 2.0, ny / 2.0]
    w.wcs.cd = [[-1e-5, 0.0], [0.0, 1e-5]]
    w.pixel_shape = (nx, ny)
    return w


def _plain_header(w):
    """A plain image-HDU header: ``NAXIS1``/``NAXIS2`` are the image."""
    nx, ny = w.pixel_shape
    header = w.to_header()
    header["NAXIS"] = 2
    header["NAXIS1"] = nx
    header["NAXIS2"] = ny
    return header


def _compressed_header(w):
    """A header shaped like an fpack tile-compressed HDU.

    ``NAXIS1``/``NAXIS2`` describe the compressed binary table (byte width, row
    count); the true image dimensions live in ``ZNAXIS1``/``ZNAXIS2``.
    """
    nx, ny = w.pixel_shape
    header = w.to_header()
    header["NAXIS"] = 2
    header["NAXIS1"] = 8
    header["NAXIS2"] = ny
    header["ZIMAGE"] = True
    header["ZNAXIS"] = 2
    header["ZNAXIS1"] = nx
    header["ZNAXIS2"] = ny
    return header


# ---------------------------------------------------------------------------
# image-shape resolution (fpack ZNAXIS vs plain NAXIS)
# ---------------------------------------------------------------------------

def test_image_shape_prefers_znaxis_for_compressed_header():
    """A compressed HDU reports ZNAXIS dims, not the binary-table NAXIS."""
    header = _compressed_header(_tan_wcs(100.0, 20.0, nx=2080, ny=4612))

    # WCS pixel_shape would wrongly report the compressed byte width (8).
    assert wcs.WCS(header).pixel_shape == (8, 4612)
    # _image_shape recovers the true image dimensions.
    assert _image_shape(header) == (2080, 4612)


def test_image_shape_falls_back_to_naxis():
    """A plain image HDU (no ZIMAGE) uses NAXIS1/NAXIS2.

    This is the live path for the workflow: ``headers-<exp>.npy`` carries the
    decompressed header astropy hands back for a tile-compressed HDU.
    """
    header = _plain_header(_tan_wcs(100.0, 20.0, nx=2080, ny=4612))

    assert _image_shape(header) == (2080, 4612)


def test_image_shape_raises_without_dimensions():
    """A header with no image dimensions raises a clear ``ValueError``."""
    header = Header()
    header["CTYPE1"] = "RA---TAN"

    with pytest.raises(ValueError, match="image dimensions"):
        _image_shape(header)


def test_compressed_header_corners_are_full_width():
    """Corners from a compressed header span the true CCD width, not 8 px."""
    header = _compressed_header(_tan_wcs(100.0, 20.0, nx=2080, ny=4612))

    ra, dec = _ccd_corners(wcs.WCS(header), _image_shape(header))

    # RA extent must reflect ~2080 px * 1e-5 deg/px * cos(dec), not 8 px.
    assert max(ra) - min(ra) > 0.01


# ---------------------------------------------------------------------------
# per-CCD corner geometry
# ---------------------------------------------------------------------------

def test_ccd_corners_returns_four_corners_around_centre():
    """The 4 corners bracket the CCD centre and cover its full pixel extent."""
    nx, ny = 2080, 4612
    w = _tan_wcs(100.0, 20.0, nx=nx, ny=ny)

    ra, dec = _ccd_corners(w, (nx, ny))

    assert len(ra) == 4
    assert len(dec) == 4
    # The footprint straddles the CRVAL centre in both coordinates.
    assert min(ra) < 100.0 < max(ra)
    assert min(dec) < 20.0 < max(dec)
    # Extent matches the CD-scaled pixel span (with edge padding).
    npt.assert_allclose(max(dec) - min(dec), ny * 1e-5, rtol=1e-3)


# ---------------------------------------------------------------------------
# nside validation
# ---------------------------------------------------------------------------

def test_check_nside_accepts_powers_of_two():
    """Powers of two for both nside values pass validation."""
    check_nside(32, 2048)


@pytest.mark.parametrize(
    "nside_coverage, nside", [(33, 2048), (32, 100), (32, 3000), (0, 2048)]
)
def test_check_nside_rejects_non_power_of_two(nside_coverage, nside):
    """A non-power-of-two nside raises ``ValueError``."""
    with pytest.raises(ValueError):
        check_nside(nside_coverage, nside)


def test_build_map_validates_nside_before_stamping():
    """``build_map`` refuses a bad nside rather than failing deep in healsparse."""
    with pytest.raises(ValueError, match="nside"):
        build_map(["100-0"], [0.0, 0.1, 0.1, 0.0], [0.0, 0.0, 0.1, 0.1],
                  32, 3000)


# ---------------------------------------------------------------------------
# build_map — nexp contract, RA-wrap and pole guards
# ---------------------------------------------------------------------------

def test_build_map_single_footprint():
    """One footprint given as flat corner lists exercises the shape guards."""
    m = build_map(
        ["100-0"], [9.9, 10.1, 10.1, 9.9], [0.0, 0.0, 0.2, 0.2], 32, 1024
    )

    assert 0 < len(m.valid_pixels) < 200
    assert m[m.valid_pixels].max() == 1


def test_build_map_nexp_counts_overlapping_exposures():
    """Two overlapping CCDs from different exposures give value 2 in overlap."""
    # Two 0.4x0.4 deg CCDs offset by 0.2 deg in RA -> a central overlap strip.
    m = build_map(
        ["100-0", "200-0"],
        [[9.8, 10.2, 10.2, 9.8], [10.0, 10.4, 10.4, 10.0]],
        [[19.8, 19.8, 20.2, 20.2], [19.8, 19.8, 20.2, 20.2]],
        32,
        1024,
    )

    values = m[m.valid_pixels]
    # The overlap is covered by both exposures (value 2); the union edges by
    # one (value 1). Both must be present; nothing exceeds 2.
    assert values.max() == 2
    assert values.min() == 1
    assert (values == 2).sum() > 0
    assert (values == 1).sum() > 0


def test_unwrap_ra_shifts_straddling_seam():
    """A polygon straddling RA=0 is put on a common (negative) branch."""
    out = unwrap_ra([359.9, 0.1, 0.1, 359.9])
    npt.assert_allclose(sorted(out), [-0.1, -0.1, 0.1, 0.1], atol=1e-6)


def test_unwrap_ra_leaves_normal_polygon_unchanged():
    """A polygon well away from the seam is returned unchanged."""
    npt.assert_allclose(
        unwrap_ra([99.9, 100.1, 100.1, 99.9]),
        [99.9, 100.1, 100.1, 99.9],
    )


def test_build_map_invariant_under_ra_shift():
    """A seam CCD's footprint matches an identical CCD shifted +10 deg in RA.

    Comparing the seam polygon's pixel count against a reference polygon of the
    same size at the same declination fails if the seam polygon were filling
    the ~360 deg complement. This is the real RA-wrap regression guard.
    """
    dec = [20.0, 20.0, 20.2, 20.2]
    m_seam = build_map(["100-0"], [359.9, 0.1, 0.1, 359.9], dec, 32, 1024)
    m_ref = build_map(["200-0"], [9.9, 10.1, 10.1, 9.9], dec, 32, 1024)

    n_seam = len(m_seam.valid_pixels)
    n_ref = len(m_ref.valid_pixels)

    # Same-size footprints at the same declination: pixel counts agree to
    # within a few boundary pixels, and are nowhere near a hemisphere.
    assert abs(n_seam - n_ref) <= max(3, int(0.1 * n_ref))
    assert 0 < n_seam < 200


def test_build_map_pole_guard_skips_polygon(capsys):
    """A polygon with |dec| near 90 deg is skipped with a warning."""
    build_map(
        ["100-0", "200-0"],
        [[10.0, 10.2, 10.2, 10.0], [10.0, 10.2, 10.2, 10.0]],
        [[89.5, 89.5, 89.7, 89.7], [20.0, 20.0, 20.2, 20.2]],
        32,
        1024,
    )

    captured = capsys.readouterr()
    assert "skipping CCD 100-0" in captured.out
    assert "Added 1 polygons" in captured.out


# ---------------------------------------------------------------------------
# median smoothing (offline; the coverage_map rule writes the raw nexp map)
# ---------------------------------------------------------------------------

def test_median_filter_fills_a_lone_hole():
    """A single low pixel inside a covered disc is pulled up to its neighbours."""
    nside = 1024
    m = hsp.HealSparseMap.make_empty(32, nside, np.uint16)
    m[hpg.query_circle(nside, 10.0, 20.0, 0.1)] = 2

    centre = hpg.angle_to_pixel(nside, 10.0, 20.0)
    m[[centre]] = 1
    assert m[centre] == 1

    assert median_filter(m)[centre] == 2
