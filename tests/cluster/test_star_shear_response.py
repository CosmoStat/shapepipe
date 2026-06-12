"""Headline guardrail: stars must carry no net metacal shear response.

Metacal deconvolves the PSF before applying the shear-response derivative.
Done right, a point source has no shape to respond to, so the per-object
response ``R = dg/dgamma`` averages to zero. A non-zero ``<R>`` is a direct
readout of PSF leakage into the response — the failure under active debugging
(Fabian handoff, 2026-06-12). Below mag 20 it is badly broken (R ~ 1, large
R1/R2 asymmetry); above mag 20 it is ~ok, so the metric masks ``20 < mag < 26``.

Target: ``<R1>``, ``<R2>`` within ``0 +/- 0.03``.

Marked ``slow`` + ``candide``: it reads real ngmix metacal catalogs on
candide. The default evaluates an EXISTING outputs dir (cheap — seconds). A
full regeneration of those outputs (the SLURM job chain,
``tile_launcher -> job_per_tile_newversion -> job_sp_14`` through the
container) is multi-hour and is wired but NOT run by the suite; opt in with
``SHAPEPIPE_REGENERATE_STAR_GRID=1``.

Run on demand:

    SHAPEPIPE_ON_CANDIDE=1 pytest tests/cluster/test_star_shear_response.py

Emits to ``tests/_artifacts/``: an R1/R2 distribution plot, a status JSON,
and a markdown summary — the seam a later GitHub Pages step publishes from.
"""

import os

import pytest

from tests.helpers.artifacts import emit_star_response_artifacts
from tests.helpers.star_response import DEFAULT_BASE_DIR, compute_star_response


pytestmark = [pytest.mark.slow, pytest.mark.candide]


def _base_dir():
    """Outputs dir to evaluate — overridable for alternate runs / regen."""
    return os.environ.get("SHAPEPIPE_STAR_GRID_OUTPUTS", DEFAULT_BASE_DIR)


@pytest.fixture(scope="module")
def star_response(artifacts_dir):
    """Compute the metric (regenerating outputs first only if asked).

    Emits artifacts as a side effect — including on a FAIL — so the status
    surface always reflects the latest run, not only green ones.
    """
    base_dir = _base_dir()

    if os.environ.get("SHAPEPIPE_REGENERATE_STAR_GRID") == "1":
        from tests.helpers.cluster import (
            slurm_available,
            submit_star_grid_chain,
            wait_for_jobs,
        )

        if not slurm_available():
            pytest.skip("regeneration requested but SLURM tools not on PATH")
        tiles = os.environ.get("SHAPEPIPE_STAR_GRID_TILES", "").split()
        if not tiles:
            pytest.skip("set SHAPEPIPE_STAR_GRID_TILES to regenerate")
        wait_for_jobs(submit_star_grid_chain(tiles))

    summary = compute_star_response(base_dir)
    status = (
        "PASS"
        if abs(summary["R1_mean"]) < summary["tolerance"]
        and abs(summary["R2_mean"]) < summary["tolerance"]
        else "FAIL"
    )
    emit_star_response_artifacts(
        artifacts_dir,
        {**summary, "status": status},
        R1=summary["R1"],
        R2=summary["R2"],
    )
    return summary


def test_star_response_r1_consistent_with_zero(star_response):
    """``<R1>`` within tolerance of zero (no net PSF leakage into R1)."""
    assert abs(star_response["R1_mean"]) < star_response["tolerance"], (
        f"<R1> = {star_response['R1_mean']:.6f} "
        f"+/- {star_response['R1_jk_err']:.6f} "
        f"exceeds +/-{star_response['tolerance']} "
        f"(N={star_response['n_obj']}, {star_response['n_tiles']} tiles): "
        "stars carry a net metacal shear response in g1 -> PSF leakage."
    )


def test_star_response_r2_consistent_with_zero(star_response):
    """``<R2>`` within tolerance of zero (no net PSF leakage into R2)."""
    assert abs(star_response["R2_mean"]) < star_response["tolerance"], (
        f"<R2> = {star_response['R2_mean']:.6f} "
        f"+/- {star_response['R2_jk_err']:.6f} "
        f"exceeds +/-{star_response['tolerance']} "
        f"(N={star_response['n_obj']}, {star_response['n_tiles']} tiles): "
        "stars carry a net metacal shear response in g2 -> PSF leakage."
    )
