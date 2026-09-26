"""Defect interpolation (``DEFECT_FILL = interpolate``): which pixels are
interpolated, and the properties of the interpolant.

:func:`interpolable_defects` picks the defect pixels that lie in a short row
or column run with clean pixels at both ends; :func:`interpolate_defects`
fills them from the nearby clean pixels with a Clough-Tocher interpolant,
averaged over the four quarter turns of the stamp and shared by every plane
(science image and metacal noise image).
"""

import numpy as np
import numpy.testing as npt
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from shapepipe.modules.ngmix_package.defect_interpolation import (
    MAX_INTERPOLATED_RUN,
    fourfold,
    interpolable_defects,
    interpolate_defects,
)


# --- interpolable_defects ---------------------------------------------------

def _oracle(defect, max_run=MAX_INTERPOLATED_RUN):
    """Brute force: walk each defect pixel's row and column run."""
    n0, n1 = defect.shape
    out = np.zeros_like(defect)
    for i, j in zip(*np.nonzero(defect)):
        for di, dj in ((0, 1), (1, 0)):
            lo, hi = (i, j), (i, j)
            while (0 <= lo[0] - di and 0 <= lo[1] - dj
                   and defect[lo[0] - di, lo[1] - dj]):
                lo = (lo[0] - di, lo[1] - dj)
            while (hi[0] + di < n0 and hi[1] + dj < n1
                   and defect[hi[0] + di, hi[1] + dj]):
                hi = (hi[0] + di, hi[1] + dj)
            length = hi[0] - lo[0] + hi[1] - lo[1] + 1
            bounded = (lo[0] - di >= 0 and lo[1] - dj >= 0
                       and hi[0] + di < n0 and hi[1] + dj < n1)
            if bounded and length <= max_run:
                out[i, j] = True
    return out


@given(
    n=st.integers(5, 21),
    density=st.floats(0.0, 0.6),
    seed=st.integers(0, 2**31 - 1),
)
@settings(max_examples=60, deadline=None)
def test_interpolable_defects_are_the_short_bounded_runs(n, density, seed):
    """A defect pixel is interpolated exactly when its row or column run is
    at most MAX_INTERPOLATED_RUN long and has clean pixels at both ends; the
    rule commutes with quarter turns.

    Failure modes: a run touching the stamp border (edge band, corner) or a
    wide hole is interpolated from one side; a 3-px bleed is left to noise;
    a clean pixel is selected; one axis is ignored, so the rule has a
    preferred direction.
    """
    defect = np.random.RandomState(seed).uniform(size=(n, n)) < density
    out = interpolable_defects(defect)
    npt.assert_array_equal(out, _oracle(defect))
    assert not out[~defect].any()
    for k in range(1, 4):
        npt.assert_array_equal(
            interpolable_defects(np.rot90(defect, k)), np.rot90(out, k)
        )


@pytest.mark.parametrize("kind,expected", [
    ("column", True), ("bleed3", True), ("finite_bleed3", True),
    ("pixel", True), ("bleed4", False), ("blob5", False),
    ("edge_band", False), ("corner", False),
])
def test_calibrated_widths_are_interpolated_and_wider_holes_are_not(
    kind, expected,
):
    """Columns, 3-px bleeds and single pixels are interpolated; 4-px bleeds,
    blobs, edge bands and corners are not."""
    n, c = 51, 25
    defect = np.zeros((n, n), dtype=bool)
    region = {
        "column": np.s_[:, c + 8],
        "bleed3": np.s_[:, c + 8:c + 11],
        "finite_bleed3": np.s_[c - 5:c + 6, c + 8:c + 11],
        "pixel": np.s_[c, c + 8],
        "bleed4": np.s_[:, c + 8:c + 12],
        "blob5": np.s_[c - 2:c + 3, c + 8:c + 13],
        "edge_band": np.s_[:, -3:],
        "corner": np.s_[:2, :2],
    }[kind]
    defect[region] = True
    out = interpolable_defects(defect)
    assert out[defect].all() if expected else not out.any()


# --- fourfold ---------------------------------------------------------------

@given(n=st.integers(2, 20), seed=st.integers(0, 2**31 - 1))
@settings(max_examples=30, deadline=None)
def test_fourfold_is_the_quarter_turn_orbit(n, seed):
    """The union contains the mask, is invariant under quarter turns, and is
    the smallest such set (idempotent)."""
    mask = np.random.RandomState(seed).uniform(size=(n, n)) < 0.1
    out = fourfold(mask)
    assert out[mask].all()
    for k in range(1, 4):
        npt.assert_array_equal(np.rot90(out, k), out)
    npt.assert_array_equal(fourfold(out), out)
    npt.assert_array_equal(
        out, mask | np.rot90(mask) | np.rot90(mask, 2) | np.rot90(mask, 3)
    )


def test_fourfold_rejects_rectangular_stamps():
    with pytest.raises(ValueError, match="square"):
        fourfold(np.zeros((11, 12), dtype=bool))


# --- interpolate_defects ----------------------------------------------------

def _mask(n=31):
    """A column, a finite 3-px bleed and a single pixel, all bounded."""
    defect = np.zeros((n, n), dtype=bool)
    defect[:, 19] = True
    defect[4:12, 7:10] = True
    defect[22, 11] = True
    return defect


@given(seed=st.integers(0, 2**31 - 1))
@settings(max_examples=20, deadline=None)
def test_interpolation_reproduces_planes_without_reading_defects(seed):
    """Linear planes are reproduced at the interpolated pixels; clean pixels,
    and defect pixels not selected, are returned untouched; the inputs are
    not modified; and no defect value (NaN or a sentinel) is ever read.

    Failure modes: defect pixels enter the support; coordinates are
    transposed or mis-rotated; the fill smooths clean pixels.
    """
    n = 31
    rng = np.random.RandomState(seed)
    rows, cols = np.indices((n, n))
    coeff = rng.uniform(-2, 2, (2, 3))
    planes = np.array([a + b * rows + c * cols for a, b, c in coeff])
    defect = _mask(n)
    defect[:, -2:] = True  # an edge band and a blob next to the column:
    defect[13:18, 21:26] = True  # defects that are not targets
    target = interpolable_defects(defect)
    assert target[15, 19]
    assert not target[13:18, 21:26].any() and not target[:, -2:].any()
    contaminated = planes.copy()
    contaminated[:, defect] = np.nan
    saved = contaminated.copy()

    out = interpolate_defects(contaminated, defect, target)

    npt.assert_allclose(out[:, target], planes[:, target], atol=1e-5)
    npt.assert_array_equal(out[:, ~target], contaminated[:, ~target])
    npt.assert_array_equal(contaminated, saved)
    contaminated[:, defect] = 1e30
    npt.assert_array_equal(
        interpolate_defects(contaminated, defect, target)[:, target],
        out[:, target],
    )


def test_interpolation_commutes_with_quarter_turns():
    """Rotating the stamp and its mask rotates the fill. A regular grid's
    Delaunay triangulation has degenerate diagonals, so a single-orientation
    interpolant does not commute; the four-orientation average does.

    Failure mode: the average over orientations is skipped, so the fill has a
    preferred direction.
    """
    n = 31
    rows, cols = np.indices((n, n))
    image = np.exp(-((rows - 15.3) ** 2 + (cols - 14.6) ** 2) / 10.0)
    image += 0.05 * np.random.RandomState(3).normal(size=(n, n))
    planes = image[None]
    defect = _mask(n)
    target = interpolable_defects(defect)
    out = interpolate_defects(planes, defect, target)
    for k in range(1, 4):
        rotated = interpolate_defects(
            np.rot90(planes, k, axes=(1, 2)), np.rot90(defect, k),
            np.rot90(target, k),
        )
        npt.assert_allclose(
            rotated, np.rot90(out, k, axes=(1, 2)), atol=1e-12, rtol=0
        )


def test_every_plane_sees_the_same_operator():
    """The fill is one linear operator applied to every plane: filling
    a * image + b * noise gives a * fill(image) + b * fill(noise), up to the
    Clough-Tocher gradient solver's tolerance.

    Failure mode: the noise image is filled differently from the science
    image (another support, triangulation or orientation set), so metacal's
    fixnoise no longer mirrors the science image's correlated noise.
    """
    n = 31
    rng = np.random.RandomState(11)
    image, noise = rng.normal(size=(2, n, n))
    defect = _mask(n)
    target = interpolable_defects(defect)
    both = interpolate_defects(np.array([image, noise]), defect, target)
    mixed = interpolate_defects(
        np.array([2.0 * image - 3.0 * noise]), defect, target
    )
    npt.assert_allclose(mixed[0], 2.0 * both[0] - 3.0 * both[1], atol=1e-5)
    alone = interpolate_defects(noise[None], defect, target)
    npt.assert_allclose(alone[0], both[1], atol=1e-5)


def test_no_target_is_a_no_op():
    planes = np.random.RandomState(5).normal(size=(2, 15, 15))
    defect = np.zeros((15, 15), dtype=bool)
    defect[:, -3:] = True
    out = interpolate_defects(planes, defect, np.zeros_like(defect))
    npt.assert_array_equal(out, planes)


def test_degenerate_support_gives_nan_not_an_error():
    """With fewer than three non-collinear clean pixels the target is NaN,
    which the caller replaces by noise."""
    defect = np.ones((5, 5), dtype=bool)
    defect[2, 1] = defect[2, 3] = False
    target = np.zeros_like(defect)
    target[2, 2] = True
    out = interpolate_defects(np.ones((1, 5, 5)), defect, target)
    assert np.isnan(out[0, 2, 2])
