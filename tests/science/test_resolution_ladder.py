"""Resolution guardrail: m and R11 along a galaxy-size ladder.

``test_mbias.py`` pins the multiplicative bias ``m`` at one comfortable
resolution. This test sweeps it: it holds the PSF fixed (Moffat,
``fwhm = 0.55``) and the injected shear fixed (``g1 = 0.02``) and walks the
galaxy half-light radius down through a ladder of size ratios
``gal_hlr / psf_fwhm`` from 1.2 (well resolved) to 0.15 (essentially a point
source). At each rung it runs the full metacal core and records ``m`` and the
response ``R11``.

**The estimator: a paired few-seed ensemble.** Each rung's ``m`` is measured
from ``N_SEEDS`` seed *pairs*: a ``+g`` and a ``−g`` arm per seed driven with a
BIT-IDENTICAL noise realisation and sub-pixel shift (same
``RandomState(seed+1000)`` into ``make_data``; galsim ``.shear()`` draws no
rng), so the ±γ difference cancels the per-seed shape noise. The single-seed
predecessor had σ_m ≈ 5e-3 — the same order as the tolerance, an assert with no
power. Pairing at the ladder's quiet ``noise = 1e-4`` operating point collapses
σ_m to a few ×1e-4 on resolved rungs (well under ``M_TOL_RESOLVED / 3``), so the
tripwire's statistical power finally backs its tolerance. The estimator is
``metacal_sim.paired_m1``, ported from the v2 breakdown-grid harness
(``run_breakdown_grid.py``) and sharing its conventions: joint bootstrap over
seed pairs for σ_m, and a degeneracy guard that marks a rung ``degenerate``
(m/m_err null) rather than dividing by an ``Rbar`` consistent with zero.

The physics, for the teaching register: as the galaxy shrinks toward the PSF
the injected shear washes out, so the recovered ``m`` grows *negative* — from
``~+0.0004`` at ratio 1.2 through zero to ``~-0.012`` at ratio 0.15. The
response ``R11`` does something less obvious: it does **not** slide monotonically.
It dips through a shallow ~0.92 minimum at mid-resolution and then **rises** back
toward 1 as the object approaches the point-source floor (a true point source is
pathological — the round Gaussian fit collapses and R stiffens). So the JSON
carries the whole shape; the assertions pin only what a *resolved* survey galaxy
must satisfy.

The guardrail: on the **resolved** rungs (``ratio >= 0.5``) the response
correction must leave ``|m| < 5e-3`` — the same few-x-1e-3 scale as
``test_mbias`` (observed max ``0.0036``). The unresolved rungs are recorded but
**not** asserted: they document where the estimator's floor lives, not a
requirement. A response mis-scale (wrong metacal step, R not applied, a constant
factor on R) shifts ``m`` by a roughly constant offset across every rung — a +5%
R error moves each resolved ``m`` by ~0.05, tripping the tripwire loudly (the
fault-injection ``response_scale`` hook, exercised in the spot-check test).

The ladder JSON (``tests/_artifacts/resolution_ladder.json``) is the publish
seam: emitted on pass *and* fail as a side effect of the module fixture, mirroring
the star-response test. Each rung now carries ``m_err`` alongside ``m`` — the
status-page rule that no number ships without its uncertainty. Existing keys are
untouched, so the artifacts.py plotter/markdown reader stays backward-compatible.

Fast + local: marked neither ``slow`` nor ``candide``; part of the inner loop.
At ``N_SEEDS = 6`` the whole ladder (8 rungs × 6 pairs × 2 arms) runs in ~40s
in-container, ~2× the single-seed predecessor — the paired power for the price
of one metacal warmup's worth of extra fits.
"""

from pathlib import Path

import pytest

from tests.helpers.artifacts import emit_resolution_ladder_artifacts
from tests.helpers.metacal_sim import METACAL_STEP, paired_m1

_ARTIFACTS_DIR = Path(__file__).resolve().parents[1] / "_artifacts"

RATIOS = [1.2, 0.9, 0.7, 0.5, 0.35, 0.25, 0.18, 0.15]
PSF_FWHM = 0.55
INJECTED_G1 = 0.02
# Paired seeds per rung. At noise=1e-4 the pairing is so effective that N=4
# already gives σ_m ≈ M_TOL/50 on the tightest resolved rung; N=6 buys a stable
# bootstrap σ_m and headroom against one unlucky rung while keeping the ladder
# at ~2× the single-seed runtime. Bump only if a rung's σ_m approaches M_TOL/3.
N_SEEDS = 6
# Base seed per rung: distinct, deterministic, never recycled across rungs so a
# rung's paired arms + its bootstrap are reproducible and independent.
SEED_BASE = 700_000
RESOLVED_MIN = 0.5
M_TOL_RESOLVED = 5e-3
# Every resolved rung must clear this before its assert has teeth. The ladder is
# built once and this is checked in the fixture so a power regression (σ_m
# creeping toward the tolerance) fails loudly rather than silently.
SIGMA_M_POWER_MAX = M_TOL_RESOLVED / 3


def _rung(ratio, response_scale=1.0):
    """One ladder rung: paired ±g1 ensemble at ``gal_hlr = ratio * psf_fwhm``.

    Returns the response-corrected multiplicative bias ``m`` and its bootstrap
    error ``m_err``, the response ``R11`` (with its error), the paired recovered
    ``g1``, the pair count, and the degeneracy flag — the paired analogue of the
    old single-seed rung, walked down the size ladder. ``response_scale`` is the
    fault-injection hook (scales Rbar; 1.0 is the clean measurement).
    """
    out = paired_m1(INJECTED_G1, gal_hlr=ratio * PSF_FWHM, psf_fwhm=PSF_FWHM,
                    n_seeds=N_SEEDS, seed0=SEED_BASE + int(round(ratio * 1000)),
                    response_scale=response_scale)
    return dict(ratio=ratio, gal_hlr=ratio * PSF_FWHM, m=out["m"],
                m_err=out["m_err"], R11=out["R11"], R11_err=out["R11_err"],
                g1_recovered=out["g1_recovered"], n_pairs=out["n_pairs"],
                degenerate=out["degenerate"], resolved=ratio >= RESOLVED_MIN)


def _provenance():
    """Container/library provenance stamped into the JSON body."""
    import galsim
    import ngmix
    import shapepipe

    return {"ngmix_version": ngmix.__version__,
            "galsim_version": galsim.__version__,
            "shapepipe_source": str(Path(shapepipe.__file__).resolve().parent)}


def _within_tol(r):
    """Rung passes the resolved tolerance: not degenerate and ``|m| < M_TOL``.

    A degenerate resolved rung (response consistent with zero) is a failure, not
    a pass — its ``m`` is null and cannot be trusted.
    """
    return not r["degenerate"] and abs(r["m"]) < M_TOL_RESOLVED


@pytest.fixture(scope="module")
def rungs():
    """Build the ladder once, emit the artifact as a side effect (pass AND fail).

    Emitting from the fixture rather than a test keeps the JSON/PNG/MD current
    whichever assertion fails first — the star-response pattern. Each rung carries
    ``m_err`` (the status-page no-number-without-uncertainty rule).
    """
    ladder = [_rung(r) for r in RATIOS]
    resolved = [r for r in ladder if r["resolved"]]
    status = "PASS" if all(_within_tol(r) for r in resolved) else "FAIL"
    rung_payload = [
        {**r, "within_tol": bool(_within_tol(r)) if r["resolved"] else None}
        for r in ladder
    ]
    emit_resolution_ladder_artifacts(
        _ARTIFACTS_DIR,
        {"test": "resolution_ladder", "status": status,
         "injected_g1": INJECTED_G1, "metacal_step": METACAL_STEP,
         "psf_fwhm": PSF_FWHM, "n_seeds": N_SEEDS,
         "resolved_ratio_threshold": RESOLVED_MIN,
         "m_tol_resolved": M_TOL_RESOLVED,
         "sigma_m_power_max": SIGMA_M_POWER_MAX,
         "rungs": rung_payload, "provenance": _provenance()},
    )
    return ladder


def test_estimator_has_power(rungs):
    """Every resolved rung's σ_m is small enough that its assert has teeth.

    A single-seed estimator gave σ_m ≈ the tolerance — the failure mode this
    upgrade fixes. Pin the power: each resolved, non-degenerate rung must carry
    ``m_err < M_TOL_RESOLVED / 3``. If pairing ever stops working (e.g. the
    shared-noise property breaks) σ_m balloons and this trips before a
    false-green ``|m|`` assert can slip through.
    """
    weak = [
        (r["ratio"], r["m_err"]) for r in rungs
        if r["resolved"] and not r["degenerate"]
        and r["m_err"] >= SIGMA_M_POWER_MAX
    ]
    assert not weak, (
        f"resolved rungs with underpowered σ_m >= {SIGMA_M_POWER_MAX:.1e}: "
        f"{[(f'{ratio}', f'{err:.5f}') for ratio, err in weak]} -> "
        "the paired estimator lost power; |m| asserts below are not trustworthy"
    )


def test_resolved_rungs_unbiased(rungs):
    """On resolved rungs (ratio >= 0.5), ``|m|`` stays below a few x 1e-3.

    The response-corrected bias must clear ``M_TOL_RESOLVED`` at every
    resolution a survey galaxy actually occupies. A response mis-scale shifts
    every rung's ``m`` by a near-constant offset and trips this loudly. A
    resolved rung that comes back degenerate (response consistent with zero) is
    also an offender — its ``m`` cannot be trusted.
    """
    offenders = [
        (r["ratio"], "degenerate" if r["degenerate"] else f"{r['m']:+.5f}")
        for r in rungs
        if r["resolved"] and not _within_tol(r)
    ]
    assert not offenders, (
        f"resolved rungs exceed |m| < {M_TOL_RESOLVED:.0e}: "
        f"{[(f'{ratio}', m) for ratio, m in offenders]} -> "
        "response-correction / deconvolution regression across the ladder"
    )


def test_response_positive_all_rungs(rungs):
    """Every rung has a non-degenerate positive response ``R11 > 0.1``.

    Mirrors ``test_mbias``: a degenerate or negative response means the
    correction cannot be trusted at that resolution.
    """
    bad = [(r["ratio"], r["R11"]) for r in rungs if not r["R11"] > 0.1]
    assert not bad, (
        f"degenerate metacal response on rungs {bad} (R11 <= 0.1); "
        "response correction cannot be trusted there"
    )


def test_fault_injection_trips_the_tripwire():
    """A +5% response mis-scale must trip the resolved-rung assert.

    The guardrail's whole value is catching a response-correction regression: if
    R is over- or under-applied by a constant factor, every resolved rung's
    ``m`` shifts by a near-constant offset (m = g/(R·γ) − 1, so scaling R by
    1.05 drives m to ≈ −0.048). This exercises that path directly via the
    ``response_scale`` hook and confirms the offset clears ``M_TOL_RESOLVED``,
    proving the upgraded (now high-power) estimator has not gone blind to the
    fault it exists to catch. One resolved rung is enough — the offset is
    ladder-wide by construction, so a single rung keeps the check fast.
    """
    faulted = _rung(0.7, response_scale=1.05)
    assert not faulted["degenerate"], "fault-injected rung should not be degenerate"
    assert abs(faulted["m"]) >= M_TOL_RESOLVED, (
        f"a +5% response mis-scale left |m| = {abs(faulted['m']):.5f} < "
        f"{M_TOL_RESOLVED:.0e} — the ladder would MISS a response regression; "
        f"expected m ≈ −0.048"
    )
