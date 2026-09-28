"""The segmentation map UberSeg sees must call THIS object self.

``uberseg_weight`` keeps the pixels whose nearest footprint carries the
object's catalogue ``NUMBER`` and zeros the rest; it never looks at the label
under the object's position. So the one thing that has to hold, for every
object in the UNIONS catalogue, is that the relabelled map carries that
object's ``NUMBER`` on its own pixels and something else everywhere else.
These tests assert exactly that, on maps small enough to read by eye —
SExtractor is not involved.
"""

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

SCRIPTS = Path(__file__).resolve().parents[2] / "workflow" / "scripts"


def _load():
    sys.path.insert(0, str(SCRIPTS))
    try:
        spec = importlib.util.spec_from_file_location(
            "_seg_relabel", SCRIPTS / "seg_relabel.py")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    finally:
        sys.path.remove(str(SCRIPTS))
    return module


seg_relabel = _load()


def _map():
    """Two SExtractor footprints, labelled 7 and 9, on a 20x20 sky."""
    seg = np.zeros((20, 20), dtype=np.int32)
    seg[2:6, 2:6] = 7
    seg[12:18, 12:18] = 9
    return seg


def test_each_object_owns_the_footprint_it_sits_in():
    """The claim is by position; the label that comes out is the NUMBER."""
    seg = _map()
    # FITS 1-indexed centres of the two footprints, in an order that is NOT
    # the label order, so a relabelling that merely renumbered would fail.
    out, counts = seg_relabel.relabel(
        seg, number=np.array([1, 2]),
        x_image=np.array([15.0, 4.0]), y_image=np.array([15.0, 4.0]))
    assert counts["matched"] == 2
    assert set(np.unique(out[seg == 9])) == {1}
    assert set(np.unique(out[seg == 7])) == {2}
    assert np.all(out[seg == 0] == 0)


def test_an_unclaimed_footprint_becomes_a_neighbour():
    """A detection no catalogue object sits in is masked, not measured."""
    seg = _map()
    out, counts = seg_relabel.relabel(
        seg, number=np.array([1]),
        x_image=np.array([4.0]), y_image=np.array([4.0]))
    assert counts["matched"] == 1
    assert set(np.unique(out[seg == 9])) == {seg_relabel.NEIGHBOUR_LABEL}
    assert seg_relabel.NEIGHBOUR_LABEL not in np.unique(out[seg == 7])


def test_an_object_on_sky_still_gets_a_self():
    """No footprint at its position -> a disc of its own NUMBER.

    Without one, uberseg_weight would zero the whole stamp and
    Ngmix._check_central_seg_label would raise on a map lacking the label.
    """
    seg = _map()
    out, counts = seg_relabel.relabel(
        seg, number=np.array([1, 5]),
        x_image=np.array([4.0, 10.0]), y_image=np.array([4.0, 10.0]),
        fallback_radius=2)
    assert counts["unclaimed"] == 1
    # Centred on the object, and nowhere else.
    assert out[9, 9] == 5
    assert np.count_nonzero(out == 5) == np.count_nonzero(
        np.add.outer(np.arange(-2, 3) ** 2, np.arange(-2, 3) ** 2) <= 4)


def test_two_objects_in_one_footprint_both_keep_a_centre():
    """The first claimant keeps the blend; the second takes back its centre."""
    seg = _map()
    out, counts = seg_relabel.relabel(
        seg, number=np.array([1, 2]),
        x_image=np.array([14.0, 16.0]), y_image=np.array([14.0, 16.0]),
        fallback_radius=1)
    assert counts == dict(matched=1, unclaimed=0, shared=1, off_image=0,
                          shared_pixel=0)
    assert out[13, 13] == 1
    assert out[15, 15] == 2
    # The claimant still holds the bulk of the footprint.
    assert np.count_nonzero(out == 1) > np.count_nonzero(out == 2)


def test_every_object_is_self_somewhere():
    """The invariant, over a randomised map: every object on the image keeps
    pixels of its own — matched, unmatched, blended or doubled."""
    rng = np.random.default_rng(0)
    seg = np.zeros((60, 60), dtype=np.int32)
    for label in range(1, 12):
        row, col = rng.integers(0, 55, size=2)
        seg[row:row + 5, col:col + 5] = label
    number = np.arange(1, 31)
    x_image = rng.uniform(1, 60, size=30)
    y_image = rng.uniform(1, 60, size=30)
    out, counts = seg_relabel.relabel(seg, number, x_image, y_image)
    # The four claim outcomes partition the catalogue; shared_pixel is a
    # separate axis, counted on top.
    assert sum(counts[k] for k in
               ("matched", "unclaimed", "shared", "off_image")) == len(number)
    assert counts["off_image"] == 0
    for num in number:
        assert np.any(out == num), f"object {num} has no self pixels"
    # And nothing outside a footprint or a disc was invented.
    assert set(np.unique(out)) <= set(number) | {0, seg_relabel.NEIGHBOUR_LABEL}


def test_two_objects_on_one_pixel_split_the_contest():
    """A pixel has one label: the lower NUMBER keeps it, the other keeps its
    disc, and the collision is counted rather than hidden."""
    seg = _map()
    out, counts = seg_relabel.relabel(
        seg, number=np.array([4, 6]), x_image=np.array([10.0, 10.2]),
        y_image=np.array([10.0, 10.1]), fallback_radius=1)
    assert counts["shared_pixel"] == 1
    assert out[9, 9] == 4
    assert np.any(out == 6)


@pytest.mark.parametrize("x, y", [(0.4, 5.0), (5.0, 61.0)])
def test_a_position_off_the_image_is_counted_not_crashed(x, y):
    seg = _map()
    out, counts = seg_relabel.relabel(
        seg, number=np.array([1]), x_image=np.array([x]),
        y_image=np.array([y]))
    assert counts["off_image"] == 1
    assert not np.any(out == 1)
