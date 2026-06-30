"""UberSeg neighbour-masking unit tests.

Covers the ``BLEND_HANDLING = uberseg`` option added to the ngmix module:

* :func:`uberseg_weight` — the Sheldon/MEDS nearest-segment Voronoi mask.
  Geometry assertions on a synthetic two-object stamp: neighbour-side pixels
  are zeroed, the surviving central core is a *single connected* region (the
  emergent "circularisation"), and the neighbour footprint is fully removed.
* :func:`prepare_ngmix_weights` — the ``noisefill`` default is byte-for-byte
  unchanged (asserted against an independent recomputation of the legacy
  three-line noise-fill on a shared RNG), while ``uberseg`` hard-masks the
  weight (weight -> 0) and leaves the image untouched.
* The error contract when ``uberseg`` is selected without a segmentation map
  (the seg-map source is plumbing-gated upstream).
"""

import numpy as np
import numpy.testing as npt
import pytest
from scipy import ndimage

from shapepipe.modules.ngmix_package.ngmix import (
    prepare_ngmix_weights,
    uberseg_weight,
)


def two_object_seg(npix=41, sep=12, r_central=3, r_neighbour=3):
    """Stamp seg map: central object (label 1) at centre, neighbour (label 2)
    offset ``sep`` pixels along the column axis. Returns (seg, centre, neigh).
    """
    seg = np.zeros((npix, npix), dtype=np.int32)
    yy, xx = np.indices((npix, npix))
    c = npix // 2
    centre = (c, c)
    neigh = (c, c + sep)
    seg[(yy - centre[0]) ** 2 + (xx - centre[1]) ** 2 <= r_central ** 2] = 1
    seg[(yy - neigh[0]) ** 2 + (xx - neigh[1]) ** 2 <= r_neighbour ** 2] = 2
    return seg, centre, neigh


def test_uberseg_zeros_neighbour_side_keeps_centre():
    """Neighbour-side pixels lose their weight; the central pixel keeps it."""
    seg, centre, neigh = two_object_seg()
    weight = np.ones_like(seg, dtype=float)

    out = uberseg_weight(weight, seg, object_number=1)

    # Central object pixel kept; deep-neighbour pixel zeroed.
    assert out[centre] == 1.0
    assert out[neigh] == 0.0
    # Every neighbour-footprint pixel is removed from the fit.
    assert np.all(out[seg == 2] == 0.0)
    # Every central-footprint pixel survives.
    assert np.all(out[seg == 1] == 1.0)
    # The input weight is not mutated in place.
    assert np.all(weight == 1.0)


def test_uberseg_core_is_single_connected_region():
    """The surviving (kept-weight) region is one connected component — the
    emergent circular core of the nearest-segment Voronoi partition."""
    seg, _, _ = two_object_seg()
    weight = np.ones_like(seg, dtype=float)

    out = uberseg_weight(weight, seg, object_number=1)

    kept = out > 0
    _, n_components = ndimage.label(kept)
    assert n_components == 1
    # The partition splits the stamp: some pixels survive, some are masked.
    assert kept.any() and (~kept).any()


def test_uberseg_partition_is_the_perpendicular_bisector():
    """With the neighbour due-right of the centre, the mask is the far-right
    half-plane beyond the footprint bisector: the left edge survives, the
    column past the neighbour is gone."""
    seg, centre, neigh = two_object_seg(npix=41, sep=12)
    weight = np.ones_like(seg, dtype=float)

    out = uberseg_weight(weight, seg, object_number=1)

    c_row, c_col = centre
    assert out[c_row, 0] == 1.0  # far side from the neighbour: kept
    assert out[c_row, -1] == 0.0  # neighbour side edge: masked


def test_uberseg_no_neighbour_is_passthrough():
    """A stamp with only the central object (or empty seg) is unchanged."""
    npix = 21
    weight = np.random.default_rng(0).random((npix, npix)) + 0.1

    # Only the central object present.
    seg = np.zeros((npix, npix), dtype=np.int32)
    seg[8:13, 8:13] = 1
    npt.assert_array_equal(uberseg_weight(weight, seg, object_number=1), weight)

    # Wholly empty seg (no detections).
    empty = np.zeros((npix, npix), dtype=np.int32)
    npt.assert_array_equal(uberseg_weight(weight, empty, object_number=1), weight)


def test_uberseg_matches_bruteforce_nearest_segment():
    """The cKDTree result equals the O(N^2) brute-force nearest-segment rule
    (Sheldon's non-C fallback) it stands in for."""
    seg, _, _ = two_object_seg(npix=25, sep=8)
    weight = np.ones_like(seg, dtype=float)

    out = uberseg_weight(weight, seg, object_number=1)

    obj = np.argwhere(seg != 0)
    labels = seg[seg != 0]
    brute = np.ones_like(weight)
    for i in range(seg.shape[0]):
        for j in range(seg.shape[1]):
            d2 = (i - obj[:, 0]) ** 2 + (j - obj[:, 1]) ** 2
            if labels[np.argmin(d2)] != 1:
                brute[i, j] = 0.0
    npt.assert_array_equal(out, brute)


# --- prepare_ngmix_weights: default unchanged, uberseg hard-masks ----------

def _gal_flag_weight(npix=41, seed=7):
    rng = np.random.default_rng(seed)
    gal = rng.normal(0.0, 3.0, (npix, npix)) + 50.0
    weight = np.ones((npix, npix))
    flag = np.zeros((npix, npix), dtype=np.int32)
    # A handful of flagged (bad) pixels — cosmic-ray-like.
    flag[5, 5] = 1
    flag[30, 12] = 2
    return gal, flag, weight


def test_noisefill_default_is_byte_identical_to_legacy():
    """The default path reproduces the legacy three-line noise-fill exactly
    (same RNG stream): masked pixels replaced by noise, weight 1/sigma^2."""
    gal, flag, weight = _gal_flag_weight()

    gal_out, w_out, noise_out = prepare_ngmix_weights(
        gal, weight, flag, np.random.RandomState(123),
    )

    # Independent recomputation of the legacy algorithm on the same stream.
    from modopt.math.stats import sigma_mad
    rng = np.random.RandomState(123)
    mask = np.copy(weight) != 0
    mask[flag != 0] = False
    sig = sigma_mad(gal)
    w_exp = mask.astype(float) / sig ** 2
    noise_exp = rng.standard_normal(gal.shape) * sig
    noise_gal = rng.standard_normal(gal.shape) * sig
    gal_exp = np.copy(gal)
    gal_exp[~mask] = noise_gal[~mask]

    npt.assert_array_equal(gal_out, gal_exp)
    npt.assert_array_equal(w_out, w_exp)
    npt.assert_array_equal(noise_out, noise_exp)


def test_uberseg_hard_masks_weight_and_leaves_image_untouched():
    """uberseg: image returned untouched, weight zeroed on neighbour-side and
    flagged pixels, positive on the central core."""
    npix = 41
    gal, flag, weight = _gal_flag_weight(npix=npix)
    seg, centre, neigh = two_object_seg(npix=npix, sep=12)

    gal_out, w_out, _ = prepare_ngmix_weights(
        gal, weight, flag, np.random.RandomState(5),
        blend_handling="uberseg", seg=seg, object_number=1,
    )

    # Image untouched under uberseg (no noise fill).
    npt.assert_array_equal(gal_out, gal)
    # Neighbour footprint hard-masked; central centre kept.
    assert np.all(w_out[seg == 2] == 0.0)
    assert w_out[centre] > 0.0
    # Flagged bad pixels remain at weight 0 (folded into the base mask).
    assert w_out[5, 5] == 0.0


def test_uberseg_requires_seg_and_object_number():
    gal, flag, weight = _gal_flag_weight()
    with pytest.raises(ValueError, match="requires a segmentation map"):
        prepare_ngmix_weights(
            gal, weight, flag, np.random.RandomState(0),
            blend_handling="uberseg",
        )


def test_unknown_blend_handling_raises():
    gal, flag, weight = _gal_flag_weight()
    with pytest.raises(ValueError, match="Unknown blend_handling"):
        prepare_ngmix_weights(
            gal, weight, flag, np.random.RandomState(0),
            blend_handling="bogus",
        )
