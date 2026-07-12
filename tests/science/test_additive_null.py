"""PSF-leakage additive null: round galaxy, elliptical PSF, recover c ~ 0.

The sibling of ``test_mbias.py``. That test injects a *known* shear into a round
PSF and checks the response-corrected ``g`` recovers it (multiplicative bias
``m``). This test injects *zero* shear into a **round** galaxy seen through an
**elliptical** PSF (``e1 = 0.05``) and checks the recovered shear stays at zero.
A round galaxy carries no shape of its own, so any non-zero additive bias
``c = g1_noshear / R11`` is PSF ellipticity leaking through an imperfect
deconvolution — the PSF-leakage / additive-bias null (shapepipe#763).

The machinery is the shared ``recover`` helper: the same ``make_data`` simulator
(now with ``psf_shear``), the same ``do_ngmix_metacal``, the same per-object
response ``R11 = (g1_1p - g1_1m)/(2*step)``. Only two things change from the
m-bias path — the galaxy is round (``shear=(0, 0)``) and the PSF is sheared
(``psf_shear=(0.05, 0)``). With the deconvolution wired correctly the true PSF is
fed as the model, so the deconvolution is unbiased and the recovered ``c`` sits
at the ~few x 1e-5 noise floor per seed. A regression that breaks the
deconvolution (deconvolves by a fitted model instead of the PSF image, drops it,
or shears the galaxy rather than the PSF) leaks the 0.05 PSF ellipticity straight
into ``c``, blowing past the tolerance in seconds.

Additive bias is statistical, so this asserts on a small deterministic seed
ensemble rather than a single draw: the twin's ``run_star_response.py`` finds
``|c| < 1e-3`` across the same 8-seed ensemble on develop (``star_response.json``,
finding ``star_response_null``), and the ``|c| < 1e-3`` tolerance clears the
observed null (mean ~5e-6, worst single seed ~5e-5) by >20x while still catching
any real leakage.

Fast + local: marked neither ``slow`` nor ``candide``; part of the inner loop.
"""

import numpy as np

from tests.helpers.metacal_sim import recover

SEEDS = list(range(8))  # deterministic ensemble; additive bias is statistical
PSF_E1 = 0.05  # true PSF ellipticity the deconvolution must remove
C_TOL = 1e-3  # |c| null (twin's published bound); mean ~5e-6, worst seed ~5e-5


def _c(seed):
    """Round galaxy through an e1=0.05 PSF; return additive bias c = g1/R11."""
    r = recover(seed, shear=(0.0, 0.0), psf_shear=(PSF_E1, 0.0))
    return r["g1"] / r["R11"], r["R11"]


def test_additive_bias_consistent_with_zero():
    """Recovered additive bias ``c`` is consistent with zero.

    A round galaxy through an ``e1 = 0.05`` PSF should yield ``c ~ 0`` once the
    deconvolution removes the PSF shape (the true PSF is fed as the model, so
    the deconvolution is unbiased and only the noise floor remains). Asserts the
    ensemble-mean ``|c|`` below the twin's published ``C_TOL``; a deconvolution
    regression that leaks PSF ellipticity pushes ``c`` toward 0.05 and trips it.
    """
    c = np.array([_c(seed)[0] for seed in SEEDS])
    c_mean = c.mean()
    assert abs(c_mean) < C_TOL, (
        f"additive bias c1 = {c_mean:+.6f} +/- {c.std():.6f} "
        f"(PSF e1 = {PSF_E1}, {len(c)} seeds) exceeds |c| < {C_TOL:.0e} -> "
        "PSF ellipticity is leaking through the deconvolution "
        "(PSF-leakage / additive-bias regression)"
    )


def test_response_nondegenerate():
    """Metacal response is non-degenerate — mirrors ``test_mbias``.

    Without a healthy ``R11`` the additive-bias ratio ``c = g1/R11`` cannot be
    trusted; observed ``R11 ~ 0.92``.
    """
    R11 = np.array([_c(seed)[1] for seed in SEEDS])
    assert np.all(np.isfinite(R11)) and R11.mean() > 0.1, (
        f"metacal response R11 = {R11.mean():.4f} (mean over seeds) is "
        "degenerate; additive-bias estimate cannot be trusted"
    )
