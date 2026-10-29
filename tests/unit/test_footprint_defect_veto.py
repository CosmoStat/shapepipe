"""D0 is exact overlap with the object's own, undilated ISO footprint."""

import numpy as np
import pytest

from shapepipe.modules.ngmix_package.ngmix import (
    defect_mask,
    footprint_defect_vetoes,
)

pytestmark = pytest.mark.decision("shape_measurement.footprint_defect_veto")


@pytest.mark.parametrize("pixel, rejected", [
    ((25, 25), True),       # on the footprint
    ((25, 26), False),      # adjacent, even though inside 7 px
    ((25, 43), True),       # on the footprint, beyond 13 px
    ((25, 42), False),      # adjacent to that distant footprint pixel
])
def test_exact_own_footprint_overlap(pixel, rejected):
    seg = np.zeros((51, 51), dtype=np.int32)
    seg[25, 25] = seg[25, 43] = 17
    seg[25, 26] = 23  # another object's footprint cannot reject this object
    defect = np.zeros_like(seg, dtype=bool)
    defect[pixel] = True
    assert footprint_defect_vetoes(defect, seg, 17) is rejected


@pytest.mark.parametrize("source", ["flag", "weight", "rms"])
def test_all_raw_defect_sources_count(source):
    seg = np.zeros((9, 9), dtype=np.int32)
    seg[4, 4] = 17
    weight, rms = np.ones((9, 9)), np.ones((9, 9))
    flag = np.zeros((9, 9), dtype=np.int32)
    {"flag": flag, "weight": weight, "rms": rms}[source][4, 4] = (
        8 if source == "flag" else 0
    )
    assert footprint_defect_vetoes(defect_mask(weight, flag, rms), seg, 17)


@pytest.mark.parametrize("border", [0, -1])
@pytest.mark.parametrize("hit", [True, False])
def test_border_touch_uses_overlap_within_stamp(border, hit):
    seg = np.zeros((51, 51), dtype=np.int32)
    seg[border, 10:13] = 17
    defect = np.zeros_like(seg, dtype=bool)
    defect[border, 12 if hit else 13] = True
    assert footprint_defect_vetoes(defect, seg, 17) is hit


@pytest.mark.parametrize("seg", [None, np.zeros((9, 9), dtype=np.int32),
                                 np.full((9, 9), 23, dtype=np.int32)])
def test_unavailable_footprint_rejects_even_clean_epoch(seg):
    assert footprint_defect_vetoes(np.zeros((9, 9), dtype=bool), seg, 17)


def test_misaligned_stamp_shapes_raise():
    with pytest.raises(ValueError, match="same shape"):
        footprint_defect_vetoes(np.zeros((9, 9), dtype=bool), np.ones((7, 7)), 1)
