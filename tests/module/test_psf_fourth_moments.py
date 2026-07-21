"""UNIT TESTS FOR PSF/STAR FOURTH-ORDER MOMENTS.

Cover the spin-2 fourth-moment combinations (``HSM_M4_1_*`` / ``HSM_M4_2_*``)
and galsim's spin-0 ``moments_rho4`` (``HSM_RHO4_*``) added to the PSFEx-interp
sky-coordinate shape measurement (:mod:`shapepipe.modules.psfex_interp_package`).

The measurement whitens the object by its own second-moment matrix, so the
predictions are analytic for elliptically-symmetric profiles:

* a round or elliptical **Gaussian** matches its own adaptive weight, so the
  whitened profile is an isotropic Gaussian: the spin-2 fourth moments vanish
  and ``rho4`` equals the Gaussian value 2.
* to get *non-zero* spin-2 fourth moments the profile must not be a single
  sheared circular profile (whitening circularises any such profile). We use a
  sum of two coaxial Gaussians of different ellipticity, whose isophote shape
  changes with radius.

The science guarantee is **frame invariance**: because the whitening and the
grid are built in the world frame (see ``_fourth_moments``), one sky object
rendered onto two differently oriented pixel grids gives the same fourth
moments. That is the key test below.

The WCS construction mirrors ``test_hsm_sky_coords.py``: a CD-matrix-like local
Jacobian with a realistic scale, rotation and small shear.
"""

import galsim
import numpy as np
import numpy.testing as npt
import pytest

from shapepipe.modules.psfex_interp_package.psfex_interp import (
    PSFExInterpolator,
    _fourth_moments,
)

# Column layout of psf_shapes / star_shapes.
M4_1, M4_2, RHO4 = 4, 5, 6

_STAMP = 101


def make_wcs(theta_deg, scale=0.187, g1=0.06, g2=-0.04):
    """A CD-matrix-like local WCS: scale, rotation ``theta_deg``, small shear."""
    th = np.deg2rad(theta_deg)
    c, s = np.cos(th), np.sin(th)
    mat = galsim.Shear(g1=g1, g2=g2).getMatrix() * scale @ np.array(
        [[c, -s], [s, c]]
    )
    return galsim.JacobianWCS(mat[0, 0], mat[0, 1], mat[1, 0], mat[1, 1])


def render(profile, wcs):
    """Render a world-frame profile onto the pixel grid defined by ``wcs``."""
    return profile.drawImage(
        nx=_STAMP, ny=_STAMP, wcs=wcs, method="no_pixel"
    ).array


def psf_row(profile, wcs):
    """Run the real ``_get_psfshapes`` on one stamp; return its shape row."""
    interp = object.__new__(PSFExInterpolator)
    interp.interp_PSFs = [render(profile, wcs)]
    interp._get_psfshapes([wcs])
    return interp.psf_shapes[0]


def star_row(profile, wcs, mask=None):
    """Run the real ``_get_starshapes`` on one star vignet; return its row."""
    stamp = render(profile, wcs)
    if mask is not None:
        stamp = np.where(mask, -1e30, stamp)
    interp = object.__new__(PSFExInterpolator)
    interp._get_starshapes(np.array([stamp]), [wcs])
    return interp.star_shapes[0]


# Sum of two coaxial Gaussians of different ellipticity: not elliptically
# symmetric, so its whitened spin-2 fourth moments are non-zero.
def composite(beta_deg=0.0):
    prof = galsim.Gaussian(sigma=0.4).shear(e1=0.5) + galsim.Gaussian(
        sigma=1.0
    ).shear(e1=0.1)
    if beta_deg:
        prof = prof.rotate(beta_deg * galsim.degrees)
    return prof


# ---------------------------------------------------------------------------
# Analytic checks on Gaussians.
# ---------------------------------------------------------------------------


def test_round_gaussian_moments_vanish_rho4_is_gaussian():
    """Round Gaussian: spin-2 fourth moments ~ 0, rho4 ~ 2 (Gaussian value)."""
    row = psf_row(galsim.Gaussian(sigma=0.6), make_wcs(0.0))
    npt.assert_allclose(row[M4_1], 0.0, atol=1e-4)
    npt.assert_allclose(row[M4_2], 0.0, atol=1e-4)
    npt.assert_allclose(row[RHO4], 2.0, rtol=1e-3)


@pytest.mark.parametrize(
    "e1, e2",
    [(0.3, 0.0), (0.0, 0.25), (0.3, 0.15), (-0.2, -0.3)],
)
def test_elliptical_gaussian_whitening_vanishes(e1, e2):
    """Elliptical Gaussian: whitening by its own 2nd moment circularises it, so
    the spin-2 fourth moments vanish regardless of (e1, e2)."""
    row = psf_row(galsim.Gaussian(sigma=0.6).shear(e1=e1, e2=e2), make_wcs(0.0))
    npt.assert_allclose(row[M4_1], 0.0, atol=1e-4)
    npt.assert_allclose(row[M4_2], 0.0, atol=1e-4)
    npt.assert_allclose(row[RHO4], 2.0, rtol=1e-3)


def test_composite_has_nonzero_spin2():
    """Coaxial two-Gaussian composite: non-Gaussian radial shape leaves a real
    spin-2 fourth moment. Axis-aligned, so only M4_1 is excited (M4_2 ~ 0)."""
    row = psf_row(composite(), make_wcs(0.0))
    assert abs(row[M4_1]) > 1e-2
    npt.assert_allclose(row[M4_2], 0.0, atol=1e-4)


# ---------------------------------------------------------------------------
# Frame invariance -- the science guarantee.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("theta", [28.0, 63.0, -45.0, 90.0])
def test_psf_fourth_moments_frame_invariant(theta):
    """One sky object rendered under two WCS orientations gives the same
    world-frame fourth moments (PSF path)."""
    prof = composite(beta_deg=30.0)  # rotated so both M4_1 and M4_2 are excited
    ref = psf_row(prof, make_wcs(0.0))
    rot = psf_row(prof, make_wcs(theta))

    assert abs(ref[M4_1]) > 1e-2 and abs(ref[M4_2]) > 1e-2  # non-trivial
    npt.assert_allclose(rot[M4_1], ref[M4_1], rtol=1e-5, atol=1e-6)
    npt.assert_allclose(rot[M4_2], ref[M4_2], rtol=1e-5, atol=1e-6)
    npt.assert_allclose(rot[RHO4], ref[RHO4], rtol=1e-5)


@pytest.mark.parametrize("theta", [28.0, -45.0])
def test_star_fourth_moments_frame_invariant(theta):
    """Same guarantee on the star path (which also handles a bad-pixel mask)."""
    prof = composite(beta_deg=30.0)
    ref = star_row(prof, make_wcs(0.0))
    rot = star_row(prof, make_wcs(theta))

    npt.assert_allclose(rot[M4_1], ref[M4_1], rtol=1e-5, atol=1e-6)
    npt.assert_allclose(rot[M4_2], ref[M4_2], rtol=1e-5, atol=1e-6)
    npt.assert_allclose(rot[RHO4], ref[RHO4], rtol=1e-5)


@pytest.mark.parametrize("theta", [0.0, 40.0])
def test_fourth_moments_invariant_under_world_shift(theta):
    """Off-center object: the sky-frame centroid subtraction makes the fourth
    moments invariant under a world-frame translation of the source.

    Every other fixture is drawn centered, so ``moments_centroid`` is ~0 and the
    recentering in ``_fourth_moments`` is a no-op. Real star vignets are only
    approximately centered, so this case shifts the source in the *world* frame
    (arcsec) and asserts the fourth moments are unchanged -- exercising the
    centroid transform that is otherwise dead in the suite.
    """
    prof = composite(beta_deg=30.0)
    wcs = make_wcs(theta)
    ref = psf_row(prof, wcs)
    shifted = psf_row(prof.shift(0.9, -1.3), wcs)  # world-frame shift, arcsec

    assert abs(ref[M4_1]) > 1e-2 and abs(ref[M4_2]) > 1e-2  # non-trivial
    npt.assert_allclose(shifted[M4_1], ref[M4_1], rtol=1e-4, atol=1e-6)
    npt.assert_allclose(shifted[M4_2], ref[M4_2], rtol=1e-4, atol=1e-6)
    npt.assert_allclose(shifted[RHO4], ref[RHO4], rtol=1e-4)


def test_psf_and_star_paths_agree():
    """PSF and star measurement of the same clean stamp agree."""
    prof = composite(beta_deg=30.0)
    wcs = make_wcs(17.0)
    p = psf_row(prof, wcs)
    s = star_row(prof, wcs)
    npt.assert_allclose(s[M4_1], p[M4_1], rtol=1e-4, atol=1e-6)
    npt.assert_allclose(s[M4_2], p[M4_2], rtol=1e-4, atol=1e-6)
    npt.assert_allclose(s[RHO4], p[RHO4], rtol=1e-4)


# ---------------------------------------------------------------------------
# Bad-pixel masking on the star path.
# ---------------------------------------------------------------------------


def test_masked_pixels_are_zeroed_before_fourth_moment_sum():
    """Star path: pixels carrying the ``-1e30`` bad-pixel sentinel are zeroed
    before the fourth-moment sum, so a masked vignet reproduces the clean
    result -- and the zeroing is load-bearing.

    Without the zeroing the sum runs over the raw ``-1e30`` sentinels and the
    fourth moments become garbage (the sum over an un-zeroed stamp is order
    ``-1e30``). We mask low-flux outskirt pixels: dropping them leaves the clean
    measurement essentially unchanged, while feeding the same sentinels through
    un-zeroed blows the result up by many orders of magnitude.
    """
    prof = composite(beta_deg=30.0)
    wcs = make_wcs(17.0)
    clean = star_row(prof, wcs)

    # A block of outskirt pixels, offset off-axis so it does not cancel in M4_1.
    mask = np.zeros((_STAMP, _STAMP), dtype=bool)
    c = _STAMP // 2
    mask[c + 18 : c + 22, c + 12 : c + 15] = True

    masked = star_row(prof, wcs, mask=mask)

    # Zeroing => the masked measurement reproduces the clean one.
    npt.assert_allclose(masked[M4_1], clean[M4_1], rtol=1e-5, atol=1e-6)
    npt.assert_allclose(masked[M4_2], clean[M4_2], rtol=1e-5, atol=1e-6)
    npt.assert_allclose(masked[RHO4], clean[RHO4], rtol=1e-5)

    # Contrast: the same sentinel stamp with the zeroing removed is garbage.
    sentinel = np.where(mask, -1e30, render(prof, wcs))
    moms = galsim.hsm.FindAdaptiveMom(
        galsim.Image(sentinel, wcs=wcs),
        badpix=galsim.Image(mask.astype(float)),
        strict=False,
        use_sky_coords=True,
    )
    raw = _fourth_moments(sentinel, moms, wcs)  # no zeroing applied
    assert abs(raw[0] - clean[M4_1]) > 1.0


# ---------------------------------------------------------------------------
# Failure handling.
# ---------------------------------------------------------------------------


def test_hsm_failure_fills_sentinels():
    """A flat stamp fails HSM: FLAG is set, spin-2 filled with 0, rho4 = -1."""
    interp = object.__new__(PSFExInterpolator)
    interp.interp_PSFs = [np.zeros((21, 21))]
    interp._get_psfshapes([make_wcs(0.0)])
    row = interp.psf_shapes[0]
    assert int(row[3]) == 1  # FLAG
    assert row[M4_1] == 0.0
    assert row[M4_2] == 0.0
    assert row[RHO4] == -1.0


def test_fourth_moments_failure_guard():
    """``_fourth_moments`` short-circuits on a ShapeData carrying an error."""
    failed = galsim.hsm.ShapeData(error_message="boom", moments_rho4=-1.0)
    assert _fourth_moments(np.zeros((5, 5)), failed, make_wcs(0.0)) == (
        0.0,
        0.0,
        -1.0,
    )
