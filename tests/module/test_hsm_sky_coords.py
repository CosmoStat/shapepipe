"""UNIT TESTS FOR HSM SKY-COORDINATE MOMENTS.

Pin the equivalence that motivates measuring PSF/star HSM adaptive moments
directly in world coordinates (``FindAdaptiveMom(use_sky_coords=True)``) instead
of measuring in the pixel frame and rotating afterwards with a WCS Jacobian.

The reference ``_transform_shape`` below is a verbatim copy of the shape-rotation
math that ``scripts/python/collate_star_cat.py`` used to apply (``transform_shape``,
removed in this change). The test draws an elliptical Gaussian on a stamp with a
nontrivial local WCS (rotation + shear + scale, matching a real CD matrix), then
checks that:

  pixel-frame moments  →  _transform_shape(local WCS Jacobian)

agrees with the sky-coordinate moments to numerical precision.

The local WCS here has positive determinant (no parity flip): production CFHT
runs decomposed with ``flip=False`` (the old ``transform_shape`` flip branch,
which is not a correct galsim reflection, never fired), so the faithful
equivalence to pin is the no-flip one.
"""

import galsim
import numpy as np
import numpy.testing as npt
import pytest


def _transform_shape(mom_list, jac):
    """Reference: the old pixel-frame → world shape rotation (verbatim).

    Copied from the removed ``collate_star_cat.transform_shape``. ``mom_list``
    is ``[g1, g2, sigma]`` measured in the pixel frame; ``jac`` is the local WCS.
    """
    scale, shear, theta, flip = jac.getDecomposition()

    sig_tmp = mom_list[2] * scale
    shape = galsim.Shear(g1=mom_list[0], g2=mom_list[1])
    if flip:
        shape = galsim.Shear(g1=-shape.g1, g2=shape.g2)
    shape = galsim.Shear(g=shape.g, beta=shape.beta + theta)
    shape = shear + shape

    return shape.g1, shape.g2, sig_tmp


# A CD-matrix-like local WCS: 0.187 arcsec/pix scale, ~28 deg rotation and a
# small shear, positive determinant (no parity flip).
_SCALE = 0.187
_THETA = 28.0 * np.pi / 180.0
_c, _s = np.cos(_THETA), np.sin(_THETA)
_SHEAR = galsim.Shear(g1=0.06, g2=-0.04)
_A = _SHEAR.getMatrix() * _SCALE @ np.array([[_c, -_s], [_s, _c]])
LOCAL_WCS = galsim.JacobianWCS(_A[0, 0], _A[0, 1], _A[1, 0], _A[1, 1])


@pytest.mark.parametrize(
    "e1, e2, sigma_pix",
    [
        (0.3, 0.1, 3.0),
        (-0.2, 0.25, 2.5),
        (0.15, -0.35, 4.0),
        (0.0, 0.0, 3.5),
    ],
)
def test_sky_coords_moments_match_jacobian_rotation(e1, e2, sigma_pix):
    """use_sky_coords=True reproduces pixel-frame moments + Jacobian rotation.

    The local WCS has positive determinant, so ``getDecomposition`` returns
    ``flip=False`` and the reference path is a pure scale+rotation+shear — the
    transformation ``use_sky_coords`` performs internally.
    """
    # Sanity: the reference decomposition sees no parity flip.
    assert LOCAL_WCS.getDecomposition()[3] is False

    # An elliptical Gaussian PSF, drawn on a unit-pixel-scale stamp; sigma is in
    # pixels because the pixel scale is 1.
    stamp = (
        galsim.Gaussian(sigma=sigma_pix)
        .shear(e1=e1, e2=e2)
        .drawImage(nx=51, ny=51, scale=1.0, method="no_pixel")
    )

    # Pixel-frame measurement, then the old Jacobian rotation.
    moms_pix = galsim.hsm.FindAdaptiveMom(
        galsim.Image(stamp.array), strict=False
    )
    g1_ref, g2_ref, sigma_ref = _transform_shape(
        [
            moms_pix.observed_shape.g1,
            moms_pix.observed_shape.g2,
            moms_pix.moments_sigma,
        ],
        LOCAL_WCS,
    )

    # Sky-coordinate measurement: same stamp, local WCS attached.
    moms_sky = galsim.hsm.FindAdaptiveMom(
        galsim.Image(stamp.array, wcs=LOCAL_WCS),
        strict=False,
        use_sky_coords=True,
    )

    npt.assert_allclose(moms_sky.observed_shape.g1, g1_ref, atol=1e-9, rtol=1e-6)
    npt.assert_allclose(moms_sky.observed_shape.g2, g2_ref, atol=1e-9, rtol=1e-6)
    npt.assert_allclose(moms_sky.moments_sigma, sigma_ref, rtol=1e-6)


def test_sky_coords_sigma_is_scaled_to_world_units():
    """The sky-frame size is the pixel size times the WCS scale (arcsec/pix).

    Pins the deliberate unit change: SIGMA_*_HSM leaves the measurement in world
    (arcsec) units rather than pixels once use_sky_coords is on.
    """
    stamp = (
        galsim.Gaussian(sigma=3.0)
        .drawImage(nx=51, ny=51, scale=1.0, method="no_pixel")
    )
    sigma_pix = galsim.hsm.FindAdaptiveMom(
        galsim.Image(stamp.array), strict=False
    ).moments_sigma
    sigma_sky = galsim.hsm.FindAdaptiveMom(
        galsim.Image(stamp.array, wcs=LOCAL_WCS),
        strict=False,
        use_sky_coords=True,
    ).moments_sigma

    npt.assert_allclose(sigma_sky, sigma_pix * _SCALE, rtol=1e-6)
