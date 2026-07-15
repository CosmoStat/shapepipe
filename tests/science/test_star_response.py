"""Additive-bias guardrail: round galaxy, elliptical PSF, recover c ~ 0.

The sibling of ``test_mbias.py``. That test injects a *known* shear into a round
PSF and checks the response-corrected ``g`` recovers it (multiplicative bias
``m``). This test injects *zero* shear into a round galaxy seen through an
**elliptical** PSF (``e1 = 0.05``) and checks the recovered shear stays at zero
(additive bias ``c = g_noshear / R``). A round galaxy carries no shape of its own,
so any non-zero ``c`` is PSF ellipticity leaking through an imperfect
deconvolution — the "star shear-response" / PSF-leakage null (shapepipe#763).

The machinery is identical to the m-bias path: the same ``make_data`` simulator
(now with ``psf_shear``), the same ``do_ngmix_metacal``, the same per-object
response ``R11 = (g1_1p - g1_1m)/(2*step)`` built from the shifted-image results.
Only two things change — the galaxy is round (``shear=(0, 0)``) and the PSF is
sheared (``psf_shear=(0.05, 0)``). With the deconvolution wired correctly the
PSF shape is removed and the recovered ``c`` sits at the few-x-1e-6 level; a
regression that breaks the deconvolution (e.g. deconvolves by a fitted model
instead of the PSF image, or drops it) leaks the 0.05 PSF ellipticity straight
into ``c``, blowing past the tolerance in seconds.

Why ``c``, not Fabian's ``R1 = R2 = 0 +/- 0.03``? Those are the *off-diagonal*
metacal response terms, which need the full 2x2 response matrix and a real
star grid to estimate; the self-contained single-object sim here most directly
supports the additive-bias ``c`` null — the pipeline-level statistic, and exactly
what the digital twin's ``run_star_response.py`` measures (it finds
``c = -5e-6``, reproduced here). Additive bias is statistical, so this asserts on
a small deterministic seed ensemble rather than a single draw: the in-container
ensemble gives ``c.mean = -5.1e-6``, ``c.std = 2.3e-5`` (worst single-seed ``c``
~5e-5), so the ``|c| < 1e-3`` tolerance clears the observed null by >20x while
still catching any real leakage.

Fast + local: marked neither ``slow`` nor ``candide``; part of the inner loop.
"""

import numpy as np


PSF_E1 = 0.05  # true PSF ellipticity the deconvolution must remove
METACAL_STEP = 0.01  # ngmix MetacalBootstrapper default shear step
SEEDS = range(8)  # deterministic ensemble; additive bias is statistical
C_TOL = 1e-3  # |c| null; ensemble mean ~5e-6, worst single seed ~5e-5


def _recover_c_with_response(seed):
    """Round galaxy through an elliptical PSF; return additive c and R11.

    Builds the per-object metacal response ``R11 = (g1_1p - g1_1m)/(2*step)``
    from the shifted-image results and forms the additive bias
    ``c = g1_noshear / R11``. This is the same response correction the m-bias
    test applies, with zero injected shear and a sheared PSF — so the recovered
    shear measures PSF leakage rather than an injected signal.
    """
    from shapepipe.modules.ngmix_package.ngmix import (
        Postage_stamp,
        do_ngmix_metacal,
        get_prior,
    )
    from shapepipe.testing.simulate import make_data

    rng = np.random.RandomState(seed)
    prior = get_prior(0.1857, rng)
    gals, psfs, _, weights, flags, jacobs = make_data(
        rng=np.random.RandomState(seed + 100),
        shear=(0.0, 0.0),
        psf_shear=(PSF_E1, 0.0),
        noise=1e-4,
        n_epochs=2,
        img_size=51,
    )
    stamp = Postage_stamp(bkg_sub=False, megacam_flip=False)
    stamp.gals, stamp.psfs, stamp.weights, stamp.flags, stamp.jacobs = (
        gals, psfs, weights, flags, jacobs,
    )
    res, _, _ = do_ngmix_metacal(stamp, prior, 1.0, rng)

    g1_noshear = res["noshear"]["g"][0]
    R11 = (res["1p"]["g"][0] - res["1m"]["g"][0]) / (2 * METACAL_STEP)
    return g1_noshear / R11, R11


def test_star_response_additive_null_below_e3():
    """Recovered additive bias ``c`` is consistent with zero in the ideal limit.

    A round galaxy through an ``e1 = 0.05`` PSF should yield ``c ~ 0`` once the
    deconvolution removes the PSF shape. Asserts the ensemble-mean ``|c|`` below
    ``C_TOL`` and a non-degenerate response ``R11``; a deconvolution regression
    that leaks PSF ellipticity pushes ``c`` toward 0.05 and trips the bound.
    """
    out = [_recover_c_with_response(seed) for seed in SEEDS]
    c = np.array([x[0] for x in out])
    R11 = np.array([x[1] for x in out])
    c_mean = c.mean()

    assert np.all(np.isfinite(R11)) and R11.min() > 0.1, (
        f"metacal response R11 = {R11.min():.4f} (min over seeds) is degenerate; "
        "additive-bias estimate cannot be trusted"
    )
    assert abs(c_mean) < C_TOL, (
        f"additive bias c1 = {c_mean:+.6f} +/- {c.std():.6f} "
        f"(PSF e1 = {PSF_E1}, {len(c)} seeds) "
        f"exceeds |c| < {C_TOL:.0e} -> PSF ellipticity is leaking through the "
        "deconvolution (PSF-leakage / additive-bias regression)"
    )
