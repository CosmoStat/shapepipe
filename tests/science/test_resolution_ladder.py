"""Resolution guardrail: m and R11 along a galaxy-size ladder.

``test_mbias.py`` pins the multiplicative bias ``m`` at one comfortable
resolution. This test sweeps it: it holds the PSF fixed (Moffat,
``fwhm = 0.55``) and the injected shear fixed (``g1 = 0.02``) and walks the
galaxy half-light radius down through a ladder of size ratios
``gal_hlr / psf_fwhm`` from 1.2 (well resolved) to 0.15 (essentially a point
source). At each rung it runs the full metacal core and records ``m`` and the
response ``R11``.

The physics, for the teaching register: as the galaxy shrinks toward the PSF
the injected shear washes out, so the recovered ``m`` grows *negative* — from
``~+0.004`` at ratio 1.2 through zero to ``~-0.012`` at ratio 0.15. The
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
R error moves each resolved ``m`` by ~0.05, tripping the tripwire loudly.

The ladder JSON (``tests/_artifacts/resolution_ladder.json``) is the publish
seam: emitted on pass *and* fail as a side effect of the module fixture, mirroring
the star-response test.

Fast + local: marked neither ``slow`` nor ``candide``; part of the inner loop.
"""

from pathlib import Path

import numpy as np
import pytest

from tests.helpers.artifacts import emit_resolution_ladder_artifacts
from tests.helpers.metacal_sim import METACAL_STEP, recover

_ARTIFACTS_DIR = Path(__file__).resolve().parents[1] / "_artifacts"

RATIOS = [1.2, 0.9, 0.7, 0.5, 0.35, 0.25, 0.18, 0.15]
PSF_FWHM = 0.55
INJECTED_G1 = 0.02
SEED = 7
RESOLVED_MIN = 0.5
M_TOL_RESOLVED = 5e-3


def _rung(ratio):
    """One ladder rung: inject g1=0.02 at ``gal_hlr = ratio * psf_fwhm``.

    Returns the response-corrected multiplicative bias ``m = g1/R11/g_inj - 1``,
    the response ``R11``, and the raw recovered ``g1`` — the same response
    correction as ``test_mbias``, walked down the size ladder.
    """
    out = recover(SEED, shear=(INJECTED_G1, 0.0),
                  gal_hlr=ratio * PSF_FWHM, psf_fwhm=PSF_FWHM)
    m = out["g1"] / out["R11"] / INJECTED_G1 - 1.0
    return dict(ratio=ratio, gal_hlr=ratio * PSF_FWHM, m=m,
                R11=out["R11"], g1_recovered=out["g1"],
                resolved=ratio >= RESOLVED_MIN)


def _provenance():
    """Container/library provenance stamped into the JSON body."""
    import galsim
    import ngmix
    import shapepipe

    return {"ngmix_version": ngmix.__version__,
            "galsim_version": galsim.__version__,
            "shapepipe_source": str(Path(shapepipe.__file__).resolve().parent)}


@pytest.fixture(scope="module")
def rungs():
    """Build the ladder once, emit the artifact as a side effect (pass AND fail).

    Emitting from the fixture rather than a test keeps the JSON/PNG/MD current
    whichever assertion fails first — the star-response pattern.
    """
    ladder = [_rung(r) for r in RATIOS]
    resolved = [r for r in ladder if r["resolved"]]
    status = "PASS" if all(abs(r["m"]) < M_TOL_RESOLVED for r in resolved) \
        else "FAIL"
    rung_payload = [
        {**r, "within_tol": bool(abs(r["m"]) < M_TOL_RESOLVED) if r["resolved"]
         else None}
        for r in ladder
    ]
    emit_resolution_ladder_artifacts(
        _ARTIFACTS_DIR,
        {"test": "resolution_ladder", "status": status,
         "injected_g1": INJECTED_G1, "metacal_step": METACAL_STEP,
         "psf_fwhm": PSF_FWHM, "seed": SEED,
         "resolved_ratio_threshold": RESOLVED_MIN,
         "m_tol_resolved": M_TOL_RESOLVED,
         "rungs": rung_payload, "provenance": _provenance()},
    )
    return ladder


def test_resolved_rungs_unbiased(rungs):
    """On resolved rungs (ratio >= 0.5), ``|m|`` stays below a few x 1e-3.

    The response-corrected bias must clear ``M_TOL_RESOLVED`` at every
    resolution a survey galaxy actually occupies. A response mis-scale shifts
    every rung's ``m`` by a near-constant offset and trips this loudly.
    """
    offenders = [
        (r["ratio"], r["m"]) for r in rungs
        if r["resolved"] and abs(r["m"]) >= M_TOL_RESOLVED
    ]
    assert not offenders, (
        f"resolved rungs exceed |m| < {M_TOL_RESOLVED:.0e}: "
        f"{[(f'{ratio}', f'{m:+.5f}') for ratio, m in offenders]} -> "
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
