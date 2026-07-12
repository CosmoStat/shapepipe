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
from sqlitedict import SqliteDict

from shapepipe.modules.ngmix_package.ngmix import (
    Ngmix,
    central_seg_label,
    prepare_ngmix_weights,
    seg_has_neighbour,
    uberseg_weight,
)


class _RecordingLogger:
    """Minimal logger that records the messages passed to ``info``."""

    def __init__(self):
        self.messages = []

    def info(self, msg, *_args, **_kwargs):
        self.messages.append(msg)


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


# --- central_seg_label: centre-pixel identification (#776 decision 3) -------

def test_central_seg_label_is_the_centre_pixel():
    """The central label is read from the centre pixel, not assumed to be 1."""
    seg, centre, _ = two_object_seg()
    assert central_seg_label(seg) == 1

    # Slide the neighbour footprint onto the stamp centre: the centre-pixel
    # label is now 2, proving the id comes from geometry, not "label 1".
    npix = seg.shape[0]
    c = npix // 2
    yy, xx = np.indices((npix, npix))
    shifted = np.zeros_like(seg)
    shifted[(yy - c) ** 2 + (xx - c) ** 2 <= 3 ** 2] = 2
    assert central_seg_label(shifted) == 2


def test_central_seg_label_sky_centre_raises():
    """A seg stamp whose centre pixel is sky (label 0) fails fast."""
    empty = np.zeros((21, 21), dtype=np.int32)
    with pytest.raises(ValueError, match="centre pixel is sky"):
        central_seg_label(empty)


# --- seg_has_neighbour: the per-object blend flag (#776 decision 4) ---------

def test_seg_has_neighbour():
    """True iff a non-central, non-zero label is present."""
    seg, _, _ = two_object_seg()
    assert seg_has_neighbour(seg, 1) is True

    # Only the central object present -> not blended.
    solo = np.zeros((21, 21), dtype=np.int32)
    solo[8:13, 8:13] = 1
    assert seg_has_neighbour(solo, 1) is False

    # Empty seg -> not blended.
    assert seg_has_neighbour(np.zeros((21, 21), dtype=np.int32), 1) is False

    # A lone footprint that *is* the central object -> not blended, even
    # though its label differs from a nominal "1".
    assert seg_has_neighbour(solo * 7, 7) is False


# --- uberseg dilation: additive neighbour-mask enlargement (#776 dec. 2) ----

def test_uberseg_dilation_grows_neighbour_mask():
    """dilate_neighbour>0 zeros a strict superset of the base (dilate=0) mask,
    and every base-masked pixel stays masked (additive-only).

    Geometric subtlety: for well-separated objects the base Voronoi cut (the
    perpendicular bisector) already masks the whole neighbour-side half-plane,
    so dilating the neighbour footprint bites only when its grown boundary
    crosses the bisector into the central Voronoi cell. A few iterations
    guarantee that crossing here (sep=8, r=3 -> ~2px footprint gap)."""
    seg, _, _ = two_object_seg(npix=41, sep=8)
    weight = np.ones_like(seg, dtype=float)

    out0 = uberseg_weight(weight, seg, object_number=1, dilate_neighbour=0)
    out3 = uberseg_weight(weight, seg, object_number=1, dilate_neighbour=3)

    # dilate=0 reproduces the validated no-dilation result byte-for-byte.
    npt.assert_array_equal(
        out0, uberseg_weight(weight, seg, object_number=1)
    )
    # Every pixel masked at dilate=0 is still masked at dilate=3 (additive).
    assert np.all(out3[out0 == 0.0] == 0.0)
    # And strictly more pixels are masked once the dilation crosses the
    # bisector into the central cell.
    assert (out3 == 0.0).sum() > (out0 == 0.0).sum()


def test_uberseg_dilation_zero_is_pure_sheldon():
    """dilate_neighbour=0 is bit-identical to the default (no-kwarg) call and
    to the O(N^2) brute-force nearest-segment rule."""
    seg, _, _ = two_object_seg(npix=25, sep=8)
    weight = np.ones_like(seg, dtype=float)

    out = uberseg_weight(weight, seg, object_number=1, dilate_neighbour=0)
    npt.assert_array_equal(out, uberseg_weight(weight, seg, object_number=1))

    obj = np.argwhere(seg != 0)
    labels = seg[seg != 0]
    brute = np.ones_like(weight)
    for i in range(seg.shape[0]):
        for j in range(seg.shape[1]):
            d2 = (i - obj[:, 0]) ** 2 + (j - obj[:, 1]) ** 2
            if labels[np.argmin(d2)] != 1:
                brute[i, j] = 0.0
    npt.assert_array_equal(out, brute)


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


def test_noisefill_ignores_seg_and_dilate_kwargs():
    """Under noisefill, passing seg / dilate_neighbour changes nothing: the
    result matches the plain default call on the same RNG stream."""
    gal, flag, weight = _gal_flag_weight()
    seg, _, _ = two_object_seg(npix=gal.shape[0], sep=12)

    base = prepare_ngmix_weights(
        gal, weight, flag, np.random.RandomState(321),
    )
    with_kwargs = prepare_ngmix_weights(
        gal, weight, flag, np.random.RandomState(321),
        seg=seg, object_number=1, dilate_neighbour=3,
    )
    for a, b in zip(with_kwargs, base):
        npt.assert_array_equal(a, b)


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


# --- _check_central_seg_label: obj_id authoritative, centre pixel diagnostic -

def _make_ngmix(tmp_path):
    """Minimal Ngmix instance (empty sqlite vignets) with a recording logger."""
    names = ("gal", "bkg", "psf", "weight", "flag", "headers")
    paths = [tmp_path / f"{name}.sqlite" for name in names]
    for path in paths:
        SqliteDict(str(path)).close()
    log = _RecordingLogger()
    ngmix = Ngmix(
        ["tile_cat.fits"] + [str(p) for p in paths[:5]],
        str(tmp_path), "-001-001", 30.0, 0.186, str(paths[5]), log,
    )
    return ngmix, log


def test_check_central_seg_label_matching_number_is_silent(tmp_path):
    """obj_id present and on the centre pixel: no warning, no raise."""
    ngmix, log = _make_ngmix(tmp_path)
    seg, _, _ = two_object_seg()  # central footprint has label 1
    log.messages.clear()  # drop the __init__ seed log
    ngmix._check_central_seg_label(seg, obj_id=1)
    assert log.messages == []


def test_check_central_seg_label_offset_centre_warns_but_proceeds(tmp_path):
    """A neighbour's label on the centre pixel (registration offset) warns but
    does not raise, as long as obj_id's footprint is present somewhere."""
    ngmix, log = _make_ngmix(tmp_path)
    seg, centre, _ = two_object_seg()
    # Overwrite the centre pixel with the neighbour label 2; obj_id=1 still
    # has its footprint elsewhere on the stamp.
    seg[centre] = 2
    ngmix._check_central_seg_label(seg, obj_id=1)
    assert any("registration offset" in m for m in log.messages)


def test_check_central_seg_label_missing_number_raises(tmp_path):
    """obj_id's footprint absent from the stamp entirely: fail fast."""
    ngmix, _ = _make_ngmix(tmp_path)
    seg, _, _ = two_object_seg()  # labels {0, 1, 2}; 99 is nowhere
    with pytest.raises(ValueError, match="contains that label"):
        ngmix._check_central_seg_label(seg, obj_id=99)


# --- FIX 3a: uberseg without seg_cat_path fails at construction -------------

def test_ngmix_init_uberseg_without_seg_cat_raises(tmp_path):
    """blend_handling='uberseg' with seg_cat_path=None raises at __init__."""
    names = ("gal", "bkg", "psf", "weight", "flag", "headers")
    paths = [tmp_path / f"{name}.sqlite" for name in names]
    for path in paths:
        SqliteDict(str(path)).close()
    with pytest.raises(ValueError, match="requires SEG_VIGNET_PATH"):
        Ngmix(
            ["tile_cat.fits"] + [str(p) for p in paths[:5]],
            str(tmp_path), "-001-001", 30.0, 0.186, str(paths[5]),
            _RecordingLogger(),
            blend_handling="uberseg",
            seg_cat_path=None,
        )


# --- FIX 3b: runner raises when SEG_VIGNET_PATH is set but missing ----------

class _FakeConfig:
    """Config stub: SEG_VIGNET_PATH is the only present option, and it resolves
    to ``seg_path`` (a path the test leaves nonexistent)."""

    def __init__(self, seg_path):
        self._seg_path = seg_path

    def getfloat(self, _sec, _key):
        return 30.0 if _key == "MAG_ZP" else 0.186

    def has_option(self, _sec, key):
        return key == "SEG_VIGNET_PATH"

    def getexpanded(self, _sec, _key):
        return self._seg_path


def test_runner_missing_seg_vignet_file_raises(tmp_path):
    """SEG_VIGNET_PATH configured but the resolved file is absent -> the runner
    fails fast with FileNotFoundError, before constructing Ngmix."""
    from shapepipe.modules.ngmix_runner import ngmix_runner

    input_file_list = ["tile_cat.fits"] + [f"in{i}.sqlite" for i in range(6)]
    seg_path = str(tmp_path / "does_not_exist.fits")
    with pytest.raises(FileNotFoundError, match="Segmentation vignet file"):
        ngmix_runner(
            input_file_list,
            {"output": str(tmp_path)},
            "-001-001",
            _FakeConfig(seg_path),
            "NGMIX_RUNNER",
            _RecordingLogger(),
        )
