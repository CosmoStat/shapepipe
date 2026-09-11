"""ShapePipe's compatibility layer between the image and the mccd release.

Two stand-ins in ``shapepipe_auxiliary_mccd`` keep MCCD fitting in the
ShapePipe image, and each is pinned here.

The starlet filters. ``starlet_filters`` stands in for
``mccd.mccd_utils.get_mr_filters`` when the image has no Sparse2D
``mr_transform`` (the ShapePipe image ships pysap without it, and every MCCD
fit then dies in ``MCCD._initialize``). The stand-in is only acceptable if it
is the same transform, so the defining properties of the B3-spline "a trous"
transform of a Dirac are pinned here: the first smoothing is the separable B3
kernel, each coarser smoothing uses the kernel dilated by two, the bands are
differences of successive smoothings, and all scales sum back to the Dirac.
(Checked once against Sparse2D ``mr_transform -t 2`` itself for 51x51 and
41x41 with 3 scales and 33x33 with 4: identical to the last bit.)

The numpy-2 stack. mccd 1.2.4 hands ``np.vstack`` a generator, which numpy 2
rejects; ``mccd.utils`` is given a numpy whose ``vstack`` accepts one.

Needs mccd and scipy, i.e. the container.
"""

import numpy as np
import pytest

aux = pytest.importorskip(
    "shapepipe.modules.mccd_package.shapepipe_auxiliary_mccd"
)

B3 = np.array([1, 4, 6, 4, 1]) / 16


def test_scales_sum_to_the_dirac():
    filters = aux.starlet_filters((51, 51), n_scales=3, coarse=True)
    assert filters.shape == (3, 51, 51)
    dirac = np.zeros((51, 51))
    dirac[25, 25] = 1
    np.testing.assert_allclose(filters.sum(axis=0), dirac, atol=1e-15)


def test_first_band_is_dirac_minus_b3_kernel():
    filters = aux.starlet_filters((51, 51), n_scales=3, coarse=True)
    c1 = np.zeros((51, 51))
    c1[23:28, 23:28] = np.outer(B3, B3)
    expected_w1 = -c1
    expected_w1[25, 25] += 1
    np.testing.assert_allclose(filters[0], expected_w1, atol=1e-15)


def test_coarse_scale_is_the_dilated_second_smoothing():
    filters = aux.starlet_filters((51, 51), n_scales=3, coarse=True)
    b3_holes = np.zeros(9)
    b3_holes[::2] = B3
    # 1D: B3 convolved with B3 dilated by 2; the 2D coarse scale is its
    # outer product, supported on +-6 pixels around the centre.
    c2_1d = np.convolve(B3, b3_holes)
    c2 = np.zeros((51, 51))
    c2[19:32, 19:32] = np.outer(c2_1d, c2_1d)
    np.testing.assert_allclose(filters[2], c2, atol=1e-15)


def test_shape_and_coarse_conventions_match_get_mr_filters():
    # get_mr_filters makes the shape odd (an even side loses one pixel) and
    # drops the coarse scale unless asked for it.
    assert aux.starlet_filters((50, 50)).shape == (2, 49, 49)
    assert aux.starlet_filters((51, 51), n_scales=4).shape == (3, 51, 51)
    trimmed = aux.starlet_filters((51, 51), coarse=True, trim=True)
    assert [f.shape for f in trimmed] == [(5, 5), (13, 13), (13, 13)]


def test_mccd_utils_vstack_takes_a_generator():
    rows = [np.arange(3), np.arange(3, 6)]
    stacked = aux.mccd.utils.np.vstack(row for row in rows)
    np.testing.assert_array_equal(stacked, np.vstack(rows))
    # Everything else is numpy itself.
    assert aux.mccd.utils.np.zeros is np.zeros
