"""Defect-fill symmetry guardrail: a filled bad column must not bias the shear.

Metacal deconvolves, shears and reconvolves the whole stamp image and never
reads the weight map. The noise-filled pixels of a defect are therefore a hole
in the galaxy light, and the reconvolution spreads that hole into the weighted
pixels around it. A hole along one detector axis, such as a bad column, is
anisotropic and gives an additive bias (Sheldon & Huff 2017, "Effects of
Missing Data"). ORing the defect mask with all its 90-degree rotations before
the fill makes the hole invariant under a 90-degree rotation, and so cancels
the bias (contract ``defect-mask-4fold-symmetrized``). A single rotation is
not enough: for a column off the centre it leaves the overlap of the column
and its image, on one diagonal, and a coherent c2.

Fixture: a round exponential galaxy through a round Moffat PSF (the
``make_data`` sim of the other science guardrails), one epoch, noise 1e-4
against flux 1000, and a one-pixel-wide flagged column through the stamp
centre or ``OFF_CENTRE`` pixels from it. Each seed runs twice: once as drawn
and once with the galaxy image rotated by 90 degrees, while the flag stays
fixed in the detector frame. The pair averages the sim's random sub-pixel
offset over a 90-degree orbit. When the filled set is invariant under a
90-degree rotation, the pair cancels that offset noise exactly. A filled set
that is not invariant keeps its bias, because the galaxy rotates and the
column does not.

The off-centre column lies inside the central-defect veto radius, so
production would drop that epoch (``epoch-central-defect-veto``). The fill's
symmetry is tested here on its own, through ``do_ngmix_metacal``.

Measured values are in the test report of the implementing PR; the positive
control below keeps the fixture's power in the suite.

Fast and local (about 2 s per case); part of the inner loop.
"""

import numpy as np
import pytest

from shapepipe.modules.ngmix_package import ngmix as ngmix_module
from shapepipe.modules.ngmix_package.ngmix import do_ngmix_metacal, get_prior
from shapepipe.testing.simulate import make_data
from tests.helpers.metacal_sim import METACAL_STEP, PIXEL_SCALE, build_stamp

SEEDS = range(8)
IMG_SIZE = 51
NOISE = 1e-4
OFF_CENTRE = 3
# |c| bound shared with test_additive_null. The symmetrized residual is
# ~1e-4 per pair; missing symmetrization gives |c1| = 0.17, and a single
# rotation gives |c2| = 0.006 off the centre.
C_TOL = 1e-3
# The unsymmetrized arm must clear C_TOL by this factor, so that the fixture
# stays sensitive to the failure it guards against.
POWER_FACTOR = 20


def _column_flag(offset):
    """One-pixel-wide flagged column ``offset`` pixels right of the centre."""
    flag = np.zeros((IMG_SIZE, IMG_SIZE), dtype=np.int32)
    flag[:, IMG_SIZE // 2 + offset] = 1
    return flag


def _metacal(seed, k, offset):
    """Metacal on one seed, with the galaxy image rotated by ``k`` * 90 deg.

    The PSF rotates with the galaxy, because it is part of the sky. The flag
    does not, because it is fixed to the detector.
    """
    rng = np.random.RandomState(seed)
    prior = get_prior(PIXEL_SCALE, rng)
    gals, psfs, sigmas, weights, _, jacobs = make_data(
        rng=np.random.RandomState(seed + 100),
        shear=(0.0, 0.0),
        noise=NOISE,
        n_epochs=1,
        img_size=IMG_SIZE,
    )
    gals = [np.rot90(g, k).copy() for g in gals]
    psfs = [np.rot90(p, k).copy() for p in psfs]
    stamp = build_stamp(
        (gals, psfs, sigmas, weights, [_column_flag(offset)], jacobs)
    )
    res, _, _ = do_ngmix_metacal(
        stamp, prior, 1.0, rng, centroid_source="hsm"
    )
    R11 = (res["1p"]["g"][0] - res["1m"]["g"][0]) / (2 * METACAL_STEP)
    R22 = (res["2p"]["g"][1] - res["2m"]["g"][1]) / (2 * METACAL_STEP)
    g1, g2 = res["noshear"]["g"]
    return g1, g2, R11, R22, res["noshear"]["flags"]


def _additive_bias(offset):
    """Orbit-paired additive bias ``(c1, c2, R11, R22)`` over ``SEEDS``."""
    runs = np.array(
        [_metacal(seed, k, offset) for seed in SEEDS for k in (0, 1)]
    )
    assert np.all(runs[:, 4] == 0), "a metacal fit failed on the fixture"
    g1, g2, R11, R22 = runs[:, :4].T
    return g1.mean() / R11.mean(), g2.mean() / R22.mean(), R11.mean(), R22.mean()


@pytest.mark.parametrize("offset", [0, OFF_CENTRE], ids=["centred", "off"])
def test_filled_bad_column_leaves_no_additive_bias(offset):
    """A flagged column through or near a round galaxy gives
    ``|c1|, |c2| < C_TOL``.

    Failure mode: the filled set is not invariant under a 90-degree rotation
    (no symmetrization, a single rotation, or a missing 180- or 270-degree
    image), so the reconvolved hole is anisotropic and leaks into the shear
    (defect-mask-4fold-symmetrized). A single rotation is invisible for the
    centred column and shows as c2 off the centre.
    """
    c1, c2, R11, R22 = _additive_bias(offset)
    assert R11 > 0.1 and R22 > 0.1, (
        f"degenerate metacal response R11 = {R11:.3f}, R22 = {R22:.3f}"
    )
    assert abs(c1) < C_TOL and abs(c2) < C_TOL, (
        f"filled bad column {offset} px off the centre biases the shear:"
        f" c1 = {c1:+.5f}, c2 = {c2:+.5f} (|c| < {C_TOL:.0e}) -> defect"
        " mask not 4-fold symmetrized before the fill"
        " (defect-mask-4fold-symmetrized)"
    )


def test_unsymmetrized_column_bias_exceeds_tolerance(monkeypatch):
    """Positive control: without symmetrization the same fixture fails.

    Failure mode: the fixture loses its power. For example, a larger galaxy or
    a changed sim default could shrink the column's bias below ``C_TOL``, so
    that the guardrail above passes whatever the code does.
    """
    monkeypatch.setattr(ngmix_module, "symmetrize_defects", lambda d: d)
    c1, _, _, _ = _additive_bias(0)
    assert abs(c1) > POWER_FACTOR * C_TOL, (
        f"unsymmetrized column bias c1 = {c1:+.5f} does not clear"
        f" {POWER_FACTOR} x C_TOL = {POWER_FACTOR * C_TOL:.0e}: the fixture"
        " no longer detects a missing symmetrization"
    )
