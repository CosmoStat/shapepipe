"""UNIT TESTS FOR PSFEX_INTERP."""

from types import SimpleNamespace

import galsim
from astropy.io import fits
import numpy as np
import numpy.testing as npt
import pytest
from hypothesis import assume, given, settings
from hypothesis import strategies as st

from shapepipe.pipeline import file_io
from shapepipe.modules.psfex_interp_package.psfex_interp import (
    PSFExInterpolator,
)


def _write_galaxy_catalogue(path, **columns):

    fits.HDUList(
        [
            fits.PrimaryHDU(),
            fits.BinTableHDU(),
            fits.BinTableHDU.from_columns(
                [
                    fits.Column(
                        name=name,
                        format="D",
                        array=np.asarray(values),
                    )
                    for name, values in columns.items()
                ]
            ),
        ]
    ).writeto(path)


def _elliptical_psf_stamps(shears, sigma=3.0, stamp=51, scale=1.0):
    """Draw elliptical Gaussian PSF stamps with given e-type distortions.

    Each ``(e1, e2)`` is applied as a galsim *distortion* (``shear(e1=, e2=)``)
    so the stamp is genuinely non-round; round stamps make ``e == g == 0`` and
    the e/g grammar fix becomes invisible.
    """
    return np.array(
        [
            galsim.Gaussian(sigma=sigma)
            .shear(e1=e1, e2=e2)
            .drawImage(nx=stamp, ny=stamp, scale=scale, method="no_pixel")
            .array
            for e1, e2 in shears
        ]
    )


@pytest.fixture(scope="module")
def psfex_interpolator(tmp_path_factory):
    """A shape-computing interpolator backed by a minimal galaxy catalogue.

    Module-scoped so it can be driven by the hypothesis property test without
    tripping the function-scoped-fixture health check; ``_get_psfshapes`` only
    reads ``interp_PSFs`` and overwrites ``psf_shapes``, so reuse is safe.
    """
    cat_path = tmp_path_factory.mktemp("psfex") / "galaxies.fits"
    _write_galaxy_catalogue(cat_path, X_IMAGE=[1.5, 2.5], Y_IMAGE=[3.5, 4.5])
    return PSFExInterpolator(
        dotpsf_path=None,
        galcat_path=str(cat_path),
        output_path=str(cat_path.parent),
        img_number="001",
        w_log=SimpleNamespace(info=lambda message: None),
        pos_params=["X_IMAGE", "Y_IMAGE"],
        get_shapes=True,
    )


def test_get_galaxy_positions_reads_requested_columns(tmp_path):

    cat_path = tmp_path / "galaxies.fits"
    _write_galaxy_catalogue(cat_path, X_IMAGE=[1.5, 2.5], Y_IMAGE=[3.5, 4.5])
    interp = PSFExInterpolator(
        dotpsf_path=None,
        galcat_path=str(cat_path),
        output_path=str(tmp_path),
        img_number="001",
        w_log=SimpleNamespace(info=lambda message: None),
        pos_params=["X_IMAGE", "Y_IMAGE"],
    )

    interp._get_galaxy_positions()

    assert interp.gal_pos.shape == (2, 2)
    npt.assert_allclose(interp.gal_pos, [[1.5, 3.5], [2.5, 4.5]])


def test_get_galaxy_positions_raises_for_missing_column(tmp_path):

    cat_path = tmp_path / "galaxies.fits"
    _write_galaxy_catalogue(cat_path, X_IMAGE=[1.5, 2.5])
    interp = PSFExInterpolator(
        dotpsf_path=None,
        galcat_path=str(cat_path),
        output_path=str(tmp_path),
        img_number="001",
        w_log=SimpleNamespace(info=lambda message: None),
        pos_params=["X_IMAGE", "Y_IMAGE"],
    )

    with pytest.raises(KeyError, match="Y_IMAGE"):
        interp._get_galaxy_positions()


# Distortion magnitude well inside |e| < 1, large enough that the e-type and
# g-type representations of the same shape differ by far more than fit noise
# (so the property has teeth), any orientation.
_e_complex = st.complex_numbers(
    min_magnitude=0.05, max_magnitude=0.6, allow_nan=False, allow_infinity=False
)


@settings(deadline=None, max_examples=25)
@given(e=_e_complex)
def test_get_psfshapes_stores_e_type_distortion(psfex_interpolator, e):
    """The stored HSM PSF shapes are e-type distortion, not g-type shear.

    ``_get_psfshapes`` reads ``moms.observed_shape.e1/.e2``; this pins it to the
    galsim *distortion* accessor and away from the reduced-shear ``.g1/.g2``.
    The check is exact against the same image's ``observed_shape`` (independent
    of how well HSM recovers the input), and the magnitude inequality is the
    tooth: for any non-round shape ``|e| > |g|`` strictly, so storing ``.g``
    instead (the pre-fix bug) would fail it.
    """
    psfex_interpolator.interp_PSFs = _elliptical_psf_stamps([(e.real, e.imag)])
    psfex_interpolator._get_psfshapes()

    moms = galsim.hsm.FindAdaptiveMom(
        galsim.Image(psfex_interpolator.interp_PSFs[0]), strict=False
    )
    e1, e2 = moms.observed_shape.e1, moms.observed_shape.e2
    g1, g2 = moms.observed_shape.g1, moms.observed_shape.g2

    mag_e = np.hypot(e1, e2)
    mag_g = np.hypot(g1, g2)
    assume(mag_e - mag_g > 1e-3)  # ensure a non-degenerate e-vs-g distinction

    # Stored values are exactly the e-type distortion components ...
    npt.assert_array_equal(psfex_interpolator.psf_shapes[0, 0], e1)
    npt.assert_array_equal(psfex_interpolator.psf_shapes[0, 1], e2)

    # ... whose magnitude strictly exceeds the g-type shear it would carry if
    # the .g/.e swap regressed.
    mag_stored = np.hypot(
        psfex_interpolator.psf_shapes[0, 0], psfex_interpolator.psf_shapes[0, 1]
    )
    npt.assert_allclose(mag_stored, mag_e)
    assert mag_stored > mag_g + 1e-4


def test_write_output_uses_hsm_grammar_and_routes_size_to_T(
    tmp_path, monkeypatch
):
    """``_write_output`` writes the HSM grammar names and stores T, not sigma.

    Captures the dict handed to ``save_as_fits`` and asserts the frozen
    grammar (``HSM_E1_PSF`` / ``HSM_E2_PSF`` / ``HSM_T_PSF`` / ``HSM_FLAG_PSF``,
    no legacy ``*_PSF_HSM``) and that the size column is routed through
    ``cs_util.size.sigma_to_T`` (``T = 2 sigma^2``) at the producer.
    """
    cat_path = tmp_path / "galaxies.fits"
    _write_galaxy_catalogue(cat_path, X_IMAGE=[1.5, 2.5], Y_IMAGE=[3.5, 4.5])
    interp = PSFExInterpolator(
        dotpsf_path=None,
        galcat_path=str(cat_path),
        output_path=str(tmp_path),
        img_number="001",
        w_log=SimpleNamespace(info=lambda message: None),
        pos_params=["X_IMAGE", "Y_IMAGE"],
        get_shapes=True,
    )
    interp.interp_PSFs = _elliptical_psf_stamps([(0.3, 0.1), (-0.2, 0.25)])
    interp._get_psfshapes()

    captured = {}

    def _capture(self, data=None, **kwargs):
        captured["data"] = data

    monkeypatch.setattr(file_io.FITSCatalogue, "save_as_fits", _capture)

    interp._write_output()
    data = captured["data"]

    assert {
        "VIGNET",
        "HSM_E1_PSF",
        "HSM_E2_PSF",
        "HSM_T_PSF",
        "HSM_FLAG_PSF",
    } <= set(data)
    for legacy in ("E1_PSF_HSM", "E2_PSF_HSM", "SIGMA_PSF_HSM", "FLAG_PSF_HSM"):
        assert legacy not in data

    # e-type distortion stored straight through from _get_psfshapes ...
    npt.assert_array_equal(data["HSM_E1_PSF"], interp.psf_shapes[:, 0])
    npt.assert_array_equal(data["HSM_E2_PSF"], interp.psf_shapes[:, 1])
    # ... and the size column is T = 2 sigma^2 (sigma_to_T), not raw sigma.
    npt.assert_allclose(data["HSM_T_PSF"], 2.0 * interp.psf_shapes[:, 2] ** 2)
    assert not np.allclose(data["HSM_T_PSF"], interp.psf_shapes[:, 2])
