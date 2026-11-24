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
    native_psf_stamps,
    sextractor_vignet_offsets,
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
def test_get_psfshapes_stores_g_type_shear(psfex_interpolator, e):
    """The stored HSM PSF shapes are g-type reduced shear, not e-type distortion.

    ``_get_psfshapes`` reads ``moms.observed_shape.g1/.g2``; this pins it to the
    galsim *reduced-shear* accessor and away from the distortion ``.e1/.e2``.
    The check is exact against the same image's ``observed_shape`` (independent
    of how well HSM recovers the input), and the magnitude inequality is the
    tooth: for any non-round shape ``|e| > |g|`` strictly, so storing ``.e``
    instead would fail it.
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

    # Stored values are exactly the g-type reduced-shear components ...
    npt.assert_array_equal(psfex_interpolator.psf_shapes[0, 0], g1)
    npt.assert_array_equal(psfex_interpolator.psf_shapes[0, 1], g2)

    # ... whose magnitude strictly subceeds the e-type distortion.
    mag_stored = np.hypot(
        psfex_interpolator.psf_shapes[0, 0], psfex_interpolator.psf_shapes[0, 1]
    )
    npt.assert_allclose(mag_stored, mag_g)
    assert mag_e > mag_stored + 1e-4


def test_write_output_uses_hsm_grammar_and_routes_size_to_T(
    tmp_path, monkeypatch
):
    """``_write_output`` writes the HSM grammar names and stores T, not sigma.

    Captures the dict handed to ``save_as_fits`` and asserts the frozen
    grammar (``HSM_G1_PSF`` / ``HSM_G2_PSF`` / ``HSM_T_PSF`` / ``HSM_FLAG_PSF``,
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
        "HSM_G1_PSF",
        "HSM_G2_PSF",
        "HSM_T_PSF",
        "HSM_FLAG_PSF",
    } <= set(data)
    for legacy in ("E1_PSF_HSM", "E2_PSF_HSM", "SIGMA_PSF_HSM", "FLAG_PSF_HSM"):
        assert legacy not in data

    # g-type distortion stored straight through from _get_psfshapes ...
    npt.assert_array_equal(data["HSM_G1_PSF"], interp.psf_shapes[:, 0])
    npt.assert_array_equal(data["HSM_G2_PSF"], interp.psf_shapes[:, 1])
    # ... and the size column is T = 2 sigma^2 (sigma_to_T), not raw sigma.
    npt.assert_allclose(data["HSM_T_PSF"], 2.0 * interp.psf_shapes[:, 2] ** 2)
    assert not np.allclose(data["HSM_T_PSF"], interp.psf_shapes[:, 2])


# --- Native-pixel rendering of PSFEx models (shapepipe#989) ----------------

# An undersampled effective PSF: a sheared Gaussian of sigma 1.1 px (FWHM
# 2.6 px, about 0.48 arcsec on MegaCam) convolved with the native pixel.
_EPSF = galsim.Convolve(
    galsim.Gaussian(sigma=1.1).shear(g1=0.06, g2=-0.03), galsim.Pixel(1.0)
)


def _native_epsf(offset=(0.0, 0.0), n=51):
    """The effective PSF sampled at native pixel centres: the analytic stamp."""
    return _EPSF.drawImage(
        nx=n, ny=n, scale=1.0, method="no_pixel", offset=offset
    ).array


def _model_grid(sampling, n_model):
    """A PSFEx-style model grid: point samples of the effective PSF.

    Samples are ``sampling`` native pixels apart, centred on index
    ``n_model // 2``, in native-pixel units (a unit-flux star's value per
    native pixel), so the grid sums to about ``1 / sampling**2``.
    """
    return (
        _EPSF.drawImage(
            nx=n_model, ny=n_model, scale=sampling, method="no_pixel"
        ).array
        / sampling**2
    )


def _write_dotpsf(path, basis, sampling):
    """Write a minimal degree-2 PSFEx ``.psf`` file with the given basis."""
    nbasis, ny, nx = basis.shape
    table = fits.BinTableHDU.from_columns(
        [
            fits.Column(
                name="PSF_MASK",
                format=f"{basis.size}E",
                dim=f"({nx}, {ny}, {nbasis})",
                array=basis[None].astype(np.float32),
            )
        ]
    )
    header = table.header
    header["POLNAXIS"] = 2
    header["POLNGRP"] = 1
    header["POLGRP1"] = 1
    header["POLGRP2"] = 1
    header["POLNAME1"] = "XWIN_IMAGE"
    header["POLNAME2"] = "YWIN_IMAGE"
    header["POLZERO1"], header["POLSCAL1"] = 1000.0, 2000.0
    header["POLZERO2"], header["POLSCAL2"] = 2000.0, 4000.0
    header["POLDEG1"] = 2
    header["PSF_SAMP"] = sampling
    header["ACCEPTED"] = 100
    header["CHI2"] = 1.0
    fits.HDUList([fits.PrimaryHDU(), table]).writeto(path)


def _interpolator(tmp_path):
    cat_path = tmp_path / "galaxies.fits"
    _write_galaxy_catalogue(cat_path, XWIN_IMAGE=[1.0], YWIN_IMAGE=[1.0])
    return PSFExInterpolator(
        dotpsf_path=None,
        galcat_path=str(cat_path),
        output_path=str(tmp_path),
        img_number="001",
        w_log=SimpleNamespace(info=lambda message: None),
        pos_params=["XWIN_IMAGE", "YWIN_IMAGE"],
    )


def _hsm(stamp):
    return galsim.hsm.FindAdaptiveMom(galsim.Image(stamp, scale=1.0))


@pytest.mark.decision("star_selection_psf.psf_stamp_rendering")
def test_unit_sampling_centred_returns_raw_grid(tmp_path):
    """At PSF_SAMP = 1 with centred stamps, interpsfex returns the raw grid.

    Bit for bit: the basis combination is the native stamp, untouched.
    """
    rng = np.random.default_rng(1)
    basis = rng.normal(size=(6, 51, 51)).astype(np.float32)
    dotpsf = tmp_path / "unit.psf"
    _write_dotpsf(dotpsf, basis, 1.0)
    pos = np.array([[100.0, 200.0], [1500.5, 3000.25]])

    psfs = _interpolator(tmp_path).interpsfex(str(dotpsf), pos)

    xs = (pos[:, 0] - 1000.0) / 2000.0
    ys = (pos[:, 1] - 2000.0) / 4000.0
    coeffs = np.column_stack(
        [np.ones_like(xs), xs, xs**2, ys, xs * ys, ys**2]
    )
    raw = np.array(
        [np.sum([c * b for c, b in zip(row, basis)], axis=0) for row in coeffs]
    )
    assert psfs.shape == (2, 51, 51)
    assert np.array_equal(psfs, raw)
    assert native_psf_stamps(raw, 1.0) is raw


@pytest.mark.decision("star_selection_psf.psf_stamp_rendering")
@pytest.mark.parametrize("sampling, n_model", [(0.5, 101), (0.7, 73), (1.0, 51)])
@pytest.mark.parametrize("offset", [(0.0, 0.0), (0.3, -0.45)])
def test_native_stamp_matches_analytic_rendering(
    tmp_path, sampling, n_model, offset
):
    """The native stamp of a model matches the effective PSF drawn natively.

    The model grid goes through a ``.psf`` file and ``interpsfex``; its
    native stamp keeps the native-pixel normalisation and matches the
    analytic native rendering at the same offset in HSM size, shape and
    centroid. Centred, the match is exact to 3e-5. Off-centre, Lanczos-3
    interpolation of point samples errs by up to 4e-3 in size and 2e-3 in
    shape (0.6e-2 px in centroid); a PSFEx model carries no such error,
    since PSFEx fits its grid through this same resampling.
    """
    basis = np.zeros((6, n_model, n_model))
    basis[0] = _model_grid(sampling, n_model)
    dotpsf = tmp_path / "model.psf"
    _write_dotpsf(dotpsf, basis, sampling)

    stamp = _interpolator(tmp_path).interpsfex(
        str(dotpsf), np.array([[1000.0, 2000.0]]), np.array([offset])
    )[0]
    truth = _native_epsf(offset)

    assert stamp.shape == (51, 51)
    npt.assert_allclose(stamp.sum(), truth.sum(), rtol=1e-3)
    tol = 5e-5 if (offset == (0.0, 0.0) and sampling in (0.5, 1.0)) else 5e-3
    got, want = _hsm(stamp), _hsm(truth)
    npt.assert_allclose(
        got.moments_sigma, want.moments_sigma, rtol=tol
    )
    npt.assert_allclose(
        [got.observed_shape.g1, got.observed_shape.g2],
        [want.observed_shape.g1, want.observed_shape.g2],
        atol=tol / 2,
    )
    npt.assert_allclose(
        [got.moments_centroid.x, got.moments_centroid.y],
        [want.moments_centroid.x, want.moments_centroid.y],
        atol=2 * tol,
    )


@pytest.mark.decision("star_selection_psf.psf_stamp_rendering")
@pytest.mark.parametrize("sampling, n_model", [(0.5, 101), (0.7, 73), (1.0, 51)])
def test_phase_moves_centroid_by_offset(sampling, n_model):
    """A phase offset moves the stamp centroid by that offset.

    To the Lanczos-3 interpolation error (under 1e-2 px); a sign or
    convention error would be off by up to a pixel.
    """
    model = _model_grid(sampling, n_model)[None]
    offsets = np.array([[0.0, 0.0], [0.4, -0.25], [-0.5, 0.5]])
    stamps = native_psf_stamps(
        np.repeat(model, len(offsets), axis=0), sampling, offsets
    )
    centre = _hsm(stamps[0]).moments_centroid
    for stamp, (dx, dy) in zip(stamps[1:], offsets[1:]):
        shift = _hsm(stamp).moments_centroid - centre
        npt.assert_allclose([shift.x, shift.y], [dx, dy], atol=1e-2)


def test_sextractor_vignet_offsets_round_the_barycentre():
    """The VIGNET centre pixel comes from X_IMAGE, not XWIN_IMAGE.

    X_IMAGE = 10.49 rounds (0-indexed, ``int(mx + 0.49999)``) to pixel 9;
    XWIN_IMAGE = 10.52 would round to 10. The offset of XWIN from the
    VIGNET centre is therefore 10.52 - 1 - 9 = 0.52 px.
    """
    pos = np.array([[10.52, 20.0]])
    barycentre = np.array([[10.49, 20.2]])
    npt.assert_allclose(
        sextractor_vignet_offsets(pos, barycentre), [[0.52, 0.0]]
    )
