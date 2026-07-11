"""Component / coordinate symmetry guardrail: shear estimation is isotropic.

Two things must hold for any shear estimator that is not silently broken.
*Isotropy*: g1 and g2 are the same physical quantity in a rotated frame, so the
pipeline must treat them identically — the metacal response is diagonal and
symmetric (``R11 ~ R22``, off-diagonals ~ 0), the multiplicative bias is the
same whether you inject in g1 or g2, and injected shear stays in the component
you put it in. *Handedness*: shear is spin-2, so a 180 deg rotation of the
image leaves the recovered ``(g1, g2)`` unchanged — the invariant the MegaCam
flip (``np.rot90(vign, k=2)`` for the upside-down CCDs) relies on.

This is the same ideal-limit sim ``test_mbias`` uses (tiny noise, Moffat PSF,
exponential galaxy, the true PSF fed as the model so deconvolution is exact),
run once per arm. It is cheap (four metacal runs, seconds) and it catches the
class of bug a single-component m-bias test cannot see: a g1<->g2 swap, a
transposed response matrix, a g2-only sign flip, or a MegaCam flip degraded to
a handedness-changing ``rot90(k=1)`` — each leaves the g1 arm looking fine and
blows up a g2 or cross-component number here.

Fast + local: unmarked (neither ``slow`` nor ``candide``), part of the inner
loop. No artifact — the m-bias / star-response surfaces own the published
status; this is a pure tripwire.
"""
import numpy as np

from tests.helpers.metacal_sim import recover

SEED = 42
INJECTED = 0.02  # per-component injected shear magnitude for the arms

# Diagonal-isotropy: |R11 - R22| observed ~1e-3; a transposed / single-component
# R breaks this by O(0.9).
R_DIAG_TOL = 5e-3
# Off-diagonal null: |R12|, |R21| observed <5e-4; a g1<->g2 leak sends one to
# O(0.9).
R_OFFDIAG_TOL = 3e-3
# Multiplicative-bias scale: matches test_mbias's M_TOL.
M_TOL = 5e-3
# Cross-component consistency: |m1 - m2| observed 1.4e-3.
M_CONSISTENCY_TOL = 3e-3
# Cross-component leakage of injected shear: observed ~2e-5 / 5e-5.
CROSS_TOL = 2e-3
# Flip invariance: 180 deg rotation leaves spin-2 shear unchanged (observed <1e-4).
FLIP_TOL = 2e-3


def _g1_arm():
    """Shear injected in g1 only, seed 42."""
    return recover(SEED, shear=(INJECTED, 0.0))


def _g2_arm():
    """Shear injected in g2 only, seed 42."""
    return recover(SEED, shear=(0.0, INJECTED))


def test_response_diagonal_isotropic():
    """R11 ~ R22: g1 and g2 respond identically.

    The metacal response to a g1 shift and a g2 shift must have the same
    magnitude on the diagonal — g1 and g2 are one quantity in a rotated frame.
    A transposed or single-component response matrix breaks the equality by
    order the response itself (~0.9).
    """
    a = _g1_arm()
    assert abs(a["R11"] - a["R22"]) < R_DIAG_TOL, (
        f"response diagonal anisotropic: R11 = {a['R11']:.4f}, "
        f"R22 = {a['R22']:.4f}, |diff| = {abs(a['R11'] - a['R22']):.4f} "
        f">= {R_DIAG_TOL:.0e} -> transposed / single-component R"
    )


def test_response_offdiagonal_null():
    """R12, R21 ~ 0 on both arms: no g1<->g2 cross-response.

    An isotropic estimator has a diagonal response. A g1<->g2 swap or leakage
    sends an off-diagonal element to order the diagonal (~0.9).
    """
    for name, arm in (("g1", _g1_arm()), ("g2", _g2_arm())):
        assert abs(arm["R12"]) < R_OFFDIAG_TOL and abs(arm["R21"]) < R_OFFDIAG_TOL, (
            f"{name} arm off-diagonal response non-null: "
            f"R12 = {arm['R12']:.5f}, R21 = {arm['R21']:.5f} "
            f"(tol {R_OFFDIAG_TOL:.0e}) -> g1<->g2 swap / leakage"
        )


def test_mbias_component_consistency():
    """m1 (g1 arm) and m2 (g2 arm) are both small and agree.

    Each arm's response-corrected multiplicative bias must sit within the
    m-bias tolerance, and the two must agree with each other. A g2-only sign
    bug leaves m1 fine and blows m2 — the cross-component check is the catch a
    g1-only test misses.
    """
    a, b = _g1_arm(), _g2_arm()
    m1 = a["g1"] / a["R11"] / INJECTED - 1.0
    m2 = b["g2"] / b["R22"] / INJECTED - 1.0
    assert abs(m1) < M_TOL and abs(m2) < M_TOL, (
        f"per-component m out of tolerance: m1 = {m1:+.5f}, m2 = {m2:+.5f} "
        f"(|m| < {M_TOL:.0e})"
    )
    assert abs(m1 - m2) < M_CONSISTENCY_TOL, (
        f"m1 and m2 inconsistent: m1 = {m1:+.5f}, m2 = {m2:+.5f}, "
        f"|m1 - m2| = {abs(m1 - m2):.5f} >= {M_CONSISTENCY_TOL:.0e} "
        f"-> a component-specific bias"
    )


def test_cross_component_recovery():
    """Injected shear stays in its own component.

    Inject in g1 -> the recovered g2 must be ~0, and vice versa. A g1<->g2 swap
    puts the shear in the wrong component and trips this.
    """
    a, b = _g1_arm(), _g2_arm()
    assert abs(a["g2"]) < CROSS_TOL, (
        f"g1-arm leaked into g2: g2 = {a['g2']:.5f} (tol {CROSS_TOL:.0e}) "
        f"-> g1<->g2 swap"
    )
    assert abs(b["g1"]) < CROSS_TOL, (
        f"g2-arm leaked into g1: g1 = {b['g1']:.5f} (tol {CROSS_TOL:.0e}) "
        f"-> g1<->g2 swap"
    )


def test_megacam_flip_consistency():
    """A 180 deg rotation of the stamp leaves recovered (g1, g2) unchanged.

    Exercises the geometry ``Ngmix.MegaCamFlip`` applies to the upside-down
    CCDs (``np.rot90(vign, k=2)``), not the flag. Shear is spin-2, so a 180 deg
    rotation must not touch it. A handedness bug (the flip degraded to
    ``rot90(k=1)``, or an axis-only flip sending g2 -> -g2) breaks the
    invariance.
    """
    ref = recover(SEED, shear=(0.02, 0.01))
    flipped = recover(SEED, shear=(0.02, 0.01), transform=lambda a: np.rot90(a, 2))
    d1, d2 = abs(ref["g1"] - flipped["g1"]), abs(ref["g2"] - flipped["g2"])
    assert d1 < FLIP_TOL and d2 < FLIP_TOL, (
        f"180 deg flip changed recovered shear: "
        f"unrotated ({ref['g1']:.5f}, {ref['g2']:.5f}) vs "
        f"rotated ({flipped['g1']:.5f}, {flipped['g2']:.5f}), "
        f"|dg1| = {d1:.5f}, |dg2| = {d2:.5f} (tol {FLIP_TOL:.0e}) "
        f"-> handedness bug (rot90 k=1 / axis-only flip)"
    )
