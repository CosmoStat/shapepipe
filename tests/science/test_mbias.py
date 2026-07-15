"""Multiplicative-bias guardrail: inject known shear, recover after R.

The headline number for any shear pipeline is the multiplicative bias ``m``,
where the recovered shear relates to the truth as ``g_rec = (1 + m) g_true``.
This test exercises the full module shape-measurement path on a controlled
simulation: a galaxy with a *known* injected shear ``g1 = 0.02`` is pushed
through ``make_ngmix_observation`` and ``do_ngmix_metacal``, the per-object
metacal response ``R11`` is built from the ``1p``/``1m`` shifted images, and
the response-corrected shear is compared to truth.

This is the ideal limit: tiny noise, Moffat PSF, exponential galaxy, the same
``make_data`` simulator the ngmix reproducibility test uses. With deconvolution
and response correction wired correctly the recovered ``g1`` lands at ~0.01996
against the injected 0.02 — ``|m| ~ 2e-3``. The assertion holds ``|m|`` below a
few x 1e-3, so a regression that breaks the response correction (e.g. drops the
deconvolution, mis-scales the metacal step) shows up as an ``m`` blowout here,
in seconds, with nothing from the cluster.

Fast + local: marked neither ``slow`` nor ``candide``; part of the inner loop.
"""

from pathlib import Path

import numpy as np
import pytest

from tests.helpers.artifacts import emit_mbias_artifacts

# The GitHub Pages publish seam (see tests/_artifacts/README.md).
_ARTIFACTS_DIR = Path(__file__).resolve().parents[1] / "_artifacts"


INJECTED_G1 = 0.02
METACAL_STEP = 0.01  # ngmix MetacalBootstrapper default shear step
M_TOL = 5e-3  # |m| in the ideal limit; recovered g1 ~ 0.01996 -> m ~ 2e-3


def _recover_g1_with_response(seed=42):
    """Inject g1=0.02, return response-corrected recovered g1.

    Builds the per-object metacal response ``R11 = (g1_1p - g1_1m)/(2*step)``
    from the shifted-image results and corrects the noshear estimate:
    ``g1_corrected = g1_noshear / R11``. This is the same response correction
    the catalog-level pipeline applies, exercised on one controlled stamp.
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
        rng=np.random.RandomState(123),
        shear=(INJECTED_G1, 0.0),
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


def test_mbias_ideal_limit_below_few_e3():
    """Response-corrected ``m`` is within a few x 1e-3 in the ideal limit.

    The controlled-shear path that already recovers g1 ~ 0.01996 vs 0.02
    injected, made into an explicit ``m`` assertion. A broken response
    correction or deconvolution pushes ``|m|`` well past the tolerance.
    """
    g1_recovered, R11 = _recover_g1_with_response()
    m = g1_recovered / INJECTED_G1 - 1.0

    # Emit the status artifact (the GitHub Pages publish seam) before asserting,
    # so the recorded status lands on pass *and* fail.
    passed = np.isfinite(R11) and R11 > 0.1 and abs(m) < M_TOL
    emit_mbias_artifacts(
        _ARTIFACTS_DIR,
        {
            "metric": "multiplicative_bias",
            "injected_g1": INJECTED_G1,
            "recovered_g1": float(g1_recovered),
            "metacal_response_R11": float(R11),
            "multiplicative_bias_m": float(m),
            "tolerance": M_TOL,
            "status": "PASS" if passed else "FAIL",
        },
    )

    assert np.isfinite(R11) and R11 > 0.1, (
        f"metacal response R11 = {R11:.4f} is degenerate; "
        "response correction cannot be trusted"
    )
    assert abs(m) < M_TOL, (
        f"multiplicative bias m = {m:.5f} "
        f"(recovered g1 = {g1_recovered:.5f} vs injected {INJECTED_G1}) "
        f"exceeds |m| < {M_TOL:.0e} in the ideal limit -> "
        "response correction / deconvolution regression"
    )
