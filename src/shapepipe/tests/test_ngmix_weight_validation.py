"""GalSim noise-injection validation of the background-RMS ngmix weights.

Empirical back pressure for the per-pixel inverse-variance weights built
from the SExtractor BACKGROUND_RMS map (``prepare_ngmix_weights`` /
``make_ngmix_observation``): galaxies with known shape are drawn with
galsim, noise of *known* per-pixel variance is injected, and the stamps
are pushed through the module's own observation/weight building and
fitter configuration. Three statistics give the test teeth:

* **chi2per**: with weights equal to the true inverse variance, the mean
  reduced chi^2 of the fit residuals is 1 to ~0.5%. This is the absolute
  check — ngmix scales its parameter covariance by reduced chi^2
  (``pcov = pcov0 * s_sq`` in ``run_leastsq``), so *pulls cannot see* a
  mis-normalised weight map, but chi2per sees every wiring break
  measured here: normalisation off x4 -> 3.99, rms/variance confusion
  -> 13.2, inverted weights -> 1.5e5, transposed map -> 13.9, scalar
  sigma_mad fallback under a gradient -> 4.8 (all vs 1.00).
* **paired shape scatter**: under a strong RMS gradient (factor 8 across
  the galaxy), inverse-variance weights must beat the old
  MegaPipe-style binary weights (``bkg_rms=None``: binarised mask /
  sigma_mad^2) on the *same* noise realisations. Measured ratio ~0.59
  direct, ~0.68 through metacal; inverted or transposed weights push it
  above 1.6.
* **pulls and mean accuracy**: reported errors are usable (pull std ~1)
  and the recovered shape is unbiased within its standard error.

The error-calibration (pull/chi2per) assertions run on direct fits of
module-built observations because metacal's fixnoise step adds a second
noise field and deconvolution-correlated noise: even with perfect
weights and flat noise the metacal noshear pull std is ~1.5 (an ngmix
fixnoise property, independent of the weight wiring), so through-metacal
pulls test ngmix, not the weights. Accuracy and the ivar-vs-binary
scatter comparison are additionally asserted end-to-end through
``do_ngmix_metacal``.
"""

import galsim
import ngmix
import numpy as np
import numpy.testing as npt
import pytest

from shapepipe.modules.ngmix_package.ngmix import (
    Postage_stamp,
    do_ngmix_metacal,
    get_prior,
    make_ngmix_observation,
    make_runners,
)

PIX_SCALE = 0.1857  # arcsec, CFIS convention
NPIX = 48
PSF_FWHM = 0.55  # arcsec
GAL_HLR = 0.5  # arcsec
GAL_FLUX = 2800.0  # S/N ~ 100 on flat noise
G_TRUE = np.array([0.03, -0.02])
SIGMA0 = 3.0  # base noise rms; deliberately != 1 so that any
# normalisation confusion (1/rms vs 1/rms^2, missing square, ...)
# shifts chi2per away from 1 instead of hiding at sigma = 1
N_REAL = 50
ACC_NSIGMA = 4.0  # accuracy tolerance in units of the standard error
# Deliberately loose: with 2N = 100 pulls the naive SE of the sample std
# is ~1/sqrt(200) ~ 0.071, but g1/g2 pulls within a realisation are
# correlated (empirical jackknife SE ~ 0.097), so 0.4 is ~4 sigma.
# Absolute error calibration is carried by the chi2per assertion, which a
# uniform 2x weight error fails (chi2per 0.50) while pulls pass (1.07).
PULL_STD_TOL = 4.0 / np.sqrt(2 * N_REAL)
CHI2_TOL = 0.05  # mean chi2per; smallest broken-wiring signature is 3.0

WCS = galsim.PixelScale(PIX_SCALE)


def rms_flat():
    """Constant noise map at SIGMA0."""
    return np.full((NPIX, NPIX), SIGMA0)


def rms_ramp():
    """Strong noise gradient: rms ramps x8 from SIGMA0 across the galaxy.

    The ramp runs through the stamp core (columns 20-28), the regime that
    motivated the rms weights: scalar noise estimates are wrong on both
    sides and flat weights average factor-64 variance differences.
    """
    x = np.broadcast_to(np.arange(NPIX, dtype=float), (NPIX, NPIX))
    return SIGMA0 * (1.0 + 7.0 * np.clip((x - 20.0) / 8.0, 0.0, 1.0))


def draw_stamp(rng, rms_map):
    """Known Gaussian galaxy x Gaussian PSF stamp with injected noise.

    Noise is drawn per pixel with the exact sigma of ``rms_map``, so the
    map fed to the module as BACKGROUND_RMS is the truth by construction.
    """
    dy, dx = rng.uniform(low=-PIX_SCALE / 2, high=PIX_SCALE / 2, size=2)
    psf = galsim.Gaussian(fwhm=PSF_FWHM)
    obj = galsim.Convolve(
        psf,
        galsim.Gaussian(half_light_radius=GAL_HLR, flux=GAL_FLUX).shear(
            g1=G_TRUE[0], g2=G_TRUE[1]
        ),
    ).shift(dx, dy)
    gal_im = obj.drawImage(nx=NPIX, ny=NPIX, wcs=WCS).array.astype(np.float64)
    psf_im = psf.drawImage(nx=NPIX, ny=NPIX, wcs=WCS).array.astype(np.float64)
    gal_im += rng.normal(size=gal_im.shape) * rms_map
    psf_im += rng.normal(scale=1e-6, size=psf_im.shape)
    return gal_im, psf_im


def fit_direct(gal, psf_im, bkg_rms, seed):
    """Fit one module-built Observation with the module's runner stack.

    Identical to ``do_ngmix_metacal`` minus the metacal image
    manipulation: same observation builder, same ``make_runners``
    configuration, plain Bootstrapper instead of MetacalBootstrapper.
    """
    rng = np.random.RandomState(seed)
    obs = make_ngmix_observation(
        gal,
        np.ones_like(gal),
        np.zeros_like(gal),
        psf_im,
        WCS.jacobian(),
        rng,
        bkg_rms=bkg_rms,
    )
    prior = get_prior(PIX_SCALE, rng)
    runner, psf_runner = make_runners(prior, GAL_FLUX, rng)
    boot = ngmix.bootstrap.Bootstrapper(
        runner=runner, psf_runner=psf_runner, ignore_failed_psf=True
    )
    obs_list = ngmix.observation.ObsList()
    obs_list.append(obs)
    return boot.go(obs_list)


def fit_metacal(gal, psf_im, bkg_rms, seed):
    """End-to-end module path: Postage_stamp -> do_ngmix_metacal noshear."""
    rng = np.random.RandomState(seed)
    stamp = Postage_stamp(bkg_sub=False, megacam_flip=False)
    stamp.gals = [gal]
    stamp.psfs = [psf_im]
    stamp.weights = [np.ones_like(gal)]
    stamp.flags = [np.zeros_like(gal)]
    stamp.jacobs = [WCS.jacobian()]
    if bkg_rms is not None:
        stamp.bkg_rms = [bkg_rms]
    res, _ = do_ngmix_metacal(stamp, get_prior(PIX_SCALE, rng), GAL_FLUX, rng)
    return res["noshear"]


def run_ensemble(fit_func, rms_map, use_ivar):
    """N_REAL seeded realisations; NaN rows mark failed fits (pairing-safe).

    Image seeds depend only on the realisation index, so ivar and binary
    runs over the same ``rms_map`` see identical noise: the scatter
    comparison is paired.
    """
    g = np.full((N_REAL, 2), np.nan)
    g_err = np.full((N_REAL, 2), np.nan)
    chi2per = np.full(N_REAL, np.nan)
    for i in range(N_REAL):
        gal, psf_im = draw_stamp(np.random.RandomState(1000 + i), rms_map)
        res = fit_func(gal, psf_im, rms_map if use_ivar else None, 5000 + i)
        if res["flags"] != 0:
            continue
        g[i] = res["g"]
        g_err[i] = np.sqrt(np.diag(res["g_cov"]))
        chi2per[i] = res.get("chi2per", np.nan)
    ok = np.isfinite(g[:, 0])
    assert ok.mean() >= 0.9, f"fit failure rate {1 - ok.mean():.0%} > 10%"
    return {
        "ok": ok,
        "resid": g - G_TRUE,
        "pulls": ((g - G_TRUE) / g_err)[ok].ravel(),
        "chi2per": chi2per[ok],
    }


@pytest.fixture(scope="module")
def direct_flat_ivar():
    return run_ensemble(fit_direct, rms_flat(), use_ivar=True)


@pytest.fixture(scope="module")
def direct_ramp_ivar():
    return run_ensemble(fit_direct, rms_ramp(), use_ivar=True)


@pytest.fixture(scope="module")
def direct_ramp_binary():
    return run_ensemble(fit_direct, rms_ramp(), use_ivar=False)


@pytest.fixture(scope="module")
def metacal_flat_ivar():
    return run_ensemble(fit_metacal, rms_flat(), use_ivar=True)


@pytest.fixture(scope="module")
def metacal_ramp_ivar():
    return run_ensemble(fit_metacal, rms_ramp(), use_ivar=True)


@pytest.fixture(scope="module")
def metacal_ramp_binary():
    return run_ensemble(fit_metacal, rms_ramp(), use_ivar=False)


def assert_unbiased(ens):
    """Mean recovered shape consistent with truth, per component.

    Tolerance is ACC_NSIGMA standard errors of the measured sample —
    nothing hand-tuned, it tightens automatically if the scatter drops.
    """
    resid = ens["resid"][ens["ok"]]
    bias = resid.mean(axis=0)
    std_err = resid.std(axis=0, ddof=1) / np.sqrt(len(resid))
    npt.assert_array_less(
        np.abs(bias),
        ACC_NSIGMA * std_err,
        err_msg=f"shape bias {bias} exceeds {ACC_NSIGMA} x SE {std_err}",
    )


def paired_scatter_ratio(ens_a, ens_b):
    """Std of shape residuals of a relative to b, on common successes.

    The std pools g1/g2 residuals about a pooled mean, so a small
    differential per-component bias inflates "scatter"; negligible at
    current magnitudes (~0.01 offsets vs ~0.07 scatter) and symmetric
    between numerator and denominator.
    """
    both = ens_a["ok"] & ens_b["ok"]
    return (
        ens_a["resid"][both].std(ddof=1) / ens_b["resid"][both].std(ddof=1)
    )


def test_direct_flat_ivar_unbiased(direct_flat_ivar):
    assert_unbiased(direct_flat_ivar)


def test_direct_ramp_ivar_unbiased(direct_ramp_ivar):
    assert_unbiased(direct_ramp_ivar)


def test_direct_ramp_binary_unbiased(direct_ramp_binary):
    """Binary weights are inefficient under a gradient, not biased."""
    assert_unbiased(direct_ramp_binary)


@pytest.mark.parametrize("case", ["direct_flat_ivar", "direct_ramp_ivar"])
def test_direct_ivar_pulls_calibrated(case, request):
    """Pull std == 1 within ~4 sigma: reported errors are real errors."""
    pull_std = request.getfixturevalue(case)["pulls"].std(ddof=1)
    assert abs(pull_std - 1.0) < PULL_STD_TOL, (
        f"{case}: pull std {pull_std:.3f} outside 1 +/- {PULL_STD_TOL:.2f}"
    )


@pytest.mark.parametrize("case", ["direct_flat_ivar", "direct_ramp_ivar"])
def test_direct_ivar_chi2_proves_absolute_normalisation(case, request):
    """Mean reduced chi^2 == 1: the weights are the true inverse variance.

    The one statistic ngmix's chi^2-rescaled covariance cannot launder.
    Every probed wiring break (normalisation, rms vs variance, inversion,
    transposition) lands at >= 3.0; correct wiring sits at 1.00 +/- 0.005.
    """
    mean_chi2 = request.getfixturevalue(case)["chi2per"].mean()
    assert abs(mean_chi2 - 1.0) < CHI2_TOL, (
        f"{case}: mean chi2per {mean_chi2:.3f} != 1 -> weight map is not"
        " the inverse of the injected per-pixel variance"
    )


def test_direct_ivar_beats_binary_under_noise_gradient(
    direct_ramp_ivar, direct_ramp_binary
):
    """The discriminating case: paired scatter ratio well below 1.

    Same noise realisations, only the weights differ. Measured ~0.59;
    binary parity gives 1.0, inverted/transposed weights give > 1.6, so
    the 0.8 threshold catches gross wiring failures (inversion, no
    downweighting). Functional-form errors that still downweight in the
    right direction — e.g. 1/rms instead of 1/rms^2 — pass the ratio
    (0.727 < 0.8) and are caught by the direct-layer chi2per assertion.
    """
    ratio = paired_scatter_ratio(direct_ramp_ivar, direct_ramp_binary)
    assert ratio < 0.8, (
        f"ivar/binary shape-scatter ratio {ratio:.3f} >= 0.8: rms weights"
        " give no precision gain under a factor-8 noise gradient"
    )


def test_metacal_flat_ivar_unbiased(metacal_flat_ivar):
    assert_unbiased(metacal_flat_ivar)


def test_metacal_ramp_ivar_unbiased(metacal_ramp_ivar):
    assert_unbiased(metacal_ramp_ivar)


def test_metacal_ivar_beats_binary_under_noise_gradient(
    metacal_ramp_ivar, metacal_ramp_binary
):
    """End-to-end through do_ngmix_metacal the ivar gain must survive.

    This pins the per-pixel fixnoise noise image: with the noise image
    drawn at a scalar median(rms) instead, fixnoise floods the low-noise
    half of the stamp and this ratio inverts to ~1.5 (ivar *worse* than
    binary). Measured ~0.68 with the per-pixel noise image.
    """
    ratio = paired_scatter_ratio(metacal_ramp_ivar, metacal_ramp_binary)
    assert ratio < 0.9, (
        f"metacal ivar/binary shape-scatter ratio {ratio:.3f} >= 0.9: the"
        " rms-weight advantage does not survive the metacal/fixnoise path"
    )
