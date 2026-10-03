"""Defect fill and the epoch cuts (ngmix module).

A defect is a stamp pixel with a nonzero flag, zero exposure weight
(off-tile pixels included) or an invalid background RMS. Before metacal,
:func:`prepare_ngmix_weights` gives every defect weight 0 and fills it the
same way whatever ``BLEND_HANDLING`` is: short defect runs take an
interpolant of the kept pixels around them, and the rest take noise at the
background RMS. Under uberseg, pixels on the neighbour side only lose their
weight, and their image values stay raw. The per-epoch cuts in
:func:`prepare_postage_stamps` act on the same defect set: the
masked-fraction cut counts it, and the central-defect veto drops an epoch
with a defect near the stamp centre, at the radius of the fill that defect
gets. The tile VIGNET's -1e30 neighbour markers are not defects: noisefill
zero-weights and noise-fills them, uberseg ignores them, and the epoch cuts
never count them.
"""

import re
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import numpy.testing as npt
import pytest
from astropy.io import fits
from astropy.wcs import WCS
from hypothesis import given, settings
from hypothesis import strategies as st
from sqlitedict import SqliteDict
from shapepipe.modules.ngmix_package import ngmix as ngmix_module

from shapepipe.modules.ngmix_package.defect_interpolation import (
    fourfold,
    interpolable_defects,
    interpolate_defects,
)
from shapepipe.modules.ngmix_package.ngmix import (
    DEFECT_WEIGHTINGS,
    EPOCH_CENTRAL_DEFECT_RADIUS,
    EPOCH_INTERPOLATED_DEFECT_RADIUS,
    EPOCH_MASKED_FRACTION_CUT,
    Ngmix,
    defect_mask,
    make_ngmix_observation,
    prepare_ngmix_weights,
    prepare_postage_stamps,
    split_tile_markers,
    uberseg_mask,
)


# --- prepare_ngmix_weights: the filled set is the defect set ---------------

@st.composite
def defect_stamps(draw):
    """Square stamp with random flagged, zero-weight and bad-RMS pixels."""
    n = draw(st.integers(min_value=5, max_value=21))
    pixels = st.lists(
        st.tuples(st.integers(0, n - 1), st.integers(0, n - 1)), max_size=n
    )
    weight = np.ones((n, n))
    flag = np.zeros((n, n), dtype=np.int32)
    bkg_rms = np.ones((n, n))
    for i, j in draw(pixels):
        weight[i, j] = 0.0
    for i, j in draw(pixels):
        flag[i, j] = draw(st.sampled_from([1, 2, 2**10]))
    for i, j in draw(pixels):
        bkg_rms[i, j] = draw(st.sampled_from([0.0, -1.0, np.nan, np.inf]))
    # Unit-noise sky with raw defect values far above it, so a filled pixel
    # is unambiguous and a leaked raw value is obvious.
    seed = draw(st.integers(0, 2**31 - 1))
    gal = np.random.RandomState(seed).normal(size=(n, n))
    valid_rms = np.isfinite(bkg_rms) & (bkg_rms > 0)
    defect = (weight == 0) | (flag != 0) | ~valid_rms
    gal[defect] = 1.0e6
    return gal, weight, flag, bkg_rms


def _uberseg_seg(n):
    """Central object on the centre pixel, neighbour footprint in a corner."""
    seg = np.zeros((n, n), dtype=np.int32)
    seg[n // 2, n // 2] = 1
    seg[:2, :2] = 2
    return seg


@given(
    stamp=defect_stamps(),
    blend_handling=st.sampled_from(["noisefill", "uberseg"]),
    seed=st.integers(0, 2**31 - 1),
)
@settings(deadline=None)
def test_filled_set_is_the_defect_set(stamp, blend_handling, seed):
    """Filled pixels are exactly the defects (flag, zero weight, bad RMS),
    and no raw defect value survives. Every other pixel keeps its raw value,
    under either BLEND_HANDLING. The weight is zero on the defects, on the
    quarter-turn copies of the interpolated ones, and under uberseg on the
    neighbour side; those copies and the neighbour side keep their light.

    Failure modes:
    * a defect source (flag, zero weight, bad RMS) is left out of the fill;
    * the filled set grows beyond the defects (for example symmetrized);
    * the zero-weight set is not the defects plus the interpolated pixels'
      quarter turns (and the uberseg neighbour side);
    * the fill is skipped under uberseg;
    * neighbour-side pixels are filled.
    """
    gal, weight, flag, bkg_rms = stamp
    n = gal.shape[0]
    kwargs = (
        dict(seg=_uberseg_seg(n), object_number=1, dilate_neighbour=1)
        if blend_handling == "uberseg"
        else {}
    )

    gal_out, w_out, _ = prepare_ngmix_weights(
        gal, weight, flag, np.random.RandomState(seed), bkg_rms=bkg_rms,
        blend_handling=blend_handling, **kwargs,
    )

    defect = (
        (weight == 0) | (flag != 0) | ~(np.isfinite(bkg_rms) & (bkg_rms > 0))
    )
    if defect.all():
        return  # fully masked: no clean pixel left to set the noise level
    filled = gal_out != gal
    npt.assert_array_equal(filled, defect, "filled set is not the defect set")
    assert np.all(np.abs(gal_out[filled]) < 100.0), "a raw defect value leaks"
    neighbour = (
        uberseg_mask(kwargs["seg"], 1, dilate_neighbour=1)
        if blend_handling == "uberseg"
        else np.zeros((n, n), dtype=bool)
    )
    # Neighbour pixels that are not defects keep their raw values
    # (filled-set equality above).
    zero = defect | neighbour | fourfold(interpolable_defects(defect))
    npt.assert_array_equal(w_out == 0.0, zero)
    npt.assert_array_equal(w_out[~zero], 1.0)


def test_uberseg_defect_in_neighbour_region_is_filled():
    """A defect pixel that also lies on the neighbour side is filled.

    Failure mode: the fill is restricted to pixels uberseg keeps, so a raw bad
    pixel on the neighbour side still reaches metacal.
    """
    n = 21
    gal = 1.0e3 + np.arange(n * n, dtype=float).reshape(n, n)
    weight = np.ones((n, n))
    flag = np.zeros((n, n), dtype=np.int32)
    flag[1, 1] = 1  # inside the neighbour footprint
    gal[1, 1] = 1.0e6
    gal_out, w_out, _ = prepare_ngmix_weights(
        gal, weight, flag, np.random.RandomState(0), bkg_rms=np.ones((n, n)),
        blend_handling="uberseg", seg=_uberseg_seg(n), object_number=1,
    )
    assert w_out[1, 1] == 0.0 and gal_out[1, 1] != gal[1, 1]
    # A neighbour-side pixel that is not a defect: zero weight, raw value.
    assert w_out[0, 3] == 0.0 and gal_out[0, 3] == gal[0, 3]


# --- prepare_postage_stamps: the fraction cut counts the defect set --------

N_STAMP = 51
RA, DEC = 150.0, 2.0


def _fake_inputs(epochs):
    """Minimal vignet / tile-catalogue stand-ins for prepare_postage_stamps.

    ``epochs`` maps ``"<exp>-<ccd>"`` to ``(flag, weight)`` or
    ``(flag, weight, bkg_rms)``. Returns ``(vignet, tile_cat, psf_obj,
    gal_obj)``.
    """
    rng = np.random.default_rng(1)
    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    wcs.wcs.crval = [RA, DEC]
    wcs.wcs.crpix = [N_STAMP / 2, N_STAMP / 2]
    wcs.wcs.cdelt = [-0.187 / 3600, 0.187 / 3600]
    header = fits.Header({"FSCALE": 1.0}).tostring()

    def per_epoch(make):
        return {k: {"VIGNET": make(k)} for k in epochs}

    with_rms = any(len(v) == 3 for v in epochs.values())
    psf_obj = per_epoch(lambda k: np.ones((N_STAMP, N_STAMP)))
    gal_obj = per_epoch(lambda k: rng.normal(0.0, 1.0, (N_STAMP, N_STAMP)))
    for epoch in gal_obj.values():
        epoch["OFFSET"] = np.zeros(2)
    vignet = SimpleNamespace(
        gal_vign_cat={"1": gal_obj},
        bkg_vign_cat=None,
        bkg_rms_vign_cat=(
            {"1": per_epoch(
                lambda k: epochs[k][2] if len(epochs[k]) == 3
                else np.ones((N_STAMP, N_STAMP))
            )}
            if with_rms
            else None
        ),
        flag_vign_cat={"1": per_epoch(lambda k: epochs[k][0])},
        weight_vign_cat={"1": per_epoch(lambda k: epochs[k][1])},
        f_wcs_file={
            k.split("-")[0]: {
                int(k.split("-")[1]): {"WCS": wcs, "header": header}
            }
            for k in epochs
        },
    )
    tile_cat = SimpleNamespace(
        vign=None, seg=None, ra=np.array([RA]), dec=np.array([DEC])
    )
    return vignet, tile_cat, psf_obj, gal_obj


def _two_sided_band(width):
    """Mask of ``width`` columns on each side of the stamp, far from the
    centre (at least 23 px for width 3)."""
    band = np.zeros((N_STAMP, N_STAMP), dtype=bool)
    band[:, :width] = True
    band[:, -width:] = True
    return band


def _epochs():
    """A clean epoch and four masked ones, all defects far from the centre.

    * band5: 5 flagged columns on one side, 9.8% (35% with its quarter turns);
    * flag6: 2 x 3 flagged columns, 11.8%;
    * dead6: the same columns at zero weight, no flags;
    * rms6: the same columns with an invalid background RMS.
    """
    clean = np.zeros((N_STAMP, N_STAMP), dtype=np.int32)
    ones = np.ones((N_STAMP, N_STAMP))
    band5 = clean.copy()
    band5[:, -5:] = 1
    two = _two_sided_band(3)
    dead = ones.copy()
    dead[two] = 0.0
    rms = ones.copy()
    rms[two] = np.nan
    return {
        "2100001-10": (clean, ones),
        "2100002-11": (band5, ones),
        "2100003-12": (two.astype(np.int32), ones),
        "2100004-13": (clean.copy(), dead),
        "2100005-14": (clean.copy(), ones, rms),
    }


def _surviving(epochs, **kwargs):
    vignet, tile_cat, psf_obj, gal_obj = _fake_inputs(epochs)
    stamp = prepare_postage_stamps(
        vignet, 1, 0, tile_cat, bkg_sub=False,
        psf_obj=psf_obj, gal_obj=gal_obj, **kwargs,
    )
    names = {id(v[0]): name for name, v in epochs.items()}
    return sorted(names[id(flag)] for flag in stamp.flags)


def test_epoch_cut_counts_the_defect_set():
    """At the 10% cut, the three 11.8% epochs are dropped, whichever defect
    source masks them, and the 9.8% band survives.

    Failure modes: the cut counts flags only (keeps dead6 and rms6), omits
    one defect source, or counts a symmetrized set (drops band5); in each
    case the cut disagrees with the set that is zero-weighted and filled
    (epoch-cut-on-defect-mask).
    """
    epochs = _epochs()
    fractions = [
        defect_mask(weight, flag, *rms).mean()
        for flag, weight, *rms in epochs.values()
    ]
    assert fractions[1] <= EPOCH_MASKED_FRACTION_CUT < min(fractions[2:])
    assert _surviving(epochs) == ["2100001-10", "2100002-11"]


# --- prepare_postage_stamps: the central-defect veto -----------------------

def _veto_epochs():
    """A clean epoch and five with one defect each, all far below the
    fraction cut: a column at EPOCH_INTERPOLATED_DEFECT_RADIUS, a column one
    pixel inside it, a single pixel at 9 px (all interpolated), and a 5-px
    bleed (noise-filled) at the interpolated radius and at
    EPOCH_CENTRAL_DEFECT_RADIUS."""
    centre = N_STAMP // 2
    ri = int(EPOCH_INTERPOLATED_DEFECT_RADIUS)
    rn = int(EPOCH_CENTRAL_DEFECT_RADIUS)
    clean = np.zeros((N_STAMP, N_STAMP), dtype=np.int32)
    ones = np.ones((N_STAMP, N_STAMP))
    column_at, column_inside, pixel, wide_at, wide_far = (
        clean.copy() for _ in range(5)
    )
    column_at[:, centre + ri] = 1
    column_inside[:, centre + ri - 1] = 1
    pixel[centre, centre + 9] = 1
    wide_at[:, centre + ri:centre + ri + 5] = 1
    wide_far[:, centre + rn:centre + rn + 5] = 1
    return {
        "2100001-10": (clean, ones),
        "2100002-11": (column_at, ones),
        "2100003-12": (column_inside, ones),
        "2100004-13": (pixel, ones),
        "2100005-14": (wide_at, ones),
        "2100006-15": (wide_far, ones),
    }


def test_the_veto_radius_follows_the_fill():
    """An interpolated defect drops the epoch inside
    EPOCH_INTERPOLATED_DEFECT_RADIUS, a noise-filled one inside
    EPOCH_CENTRAL_DEFECT_RADIUS; a defect at its radius is kept.

    Failure modes: the veto is skipped, so a filled hole in the object's
    light enters the fit; the interpolated radius is applied to noise-filled
    pixels (a wide hole next to the object survives), or not applied at all
    (a 9-px pixel is dropped); the boundary is inclusive
    (veto-radius-follows-the-fill).
    """
    assert _surviving(_veto_epochs()) == [
        "2100001-10", "2100002-11", "2100004-13", "2100006-15"
    ]


# --- prepare_postage_stamps: per-epoch OFFSET ------------------------------

def test_each_surviving_epoch_carries_its_own_offset():
    """``stamp.offsets`` holds each surviving epoch's vignette OFFSET, in the
    order of ``stamp.flags``.

    Failure mode: the offset is dropped or read from another epoch, so the
    default "wcs" centroid raises or puts the Jacobian origin off the object.
    """
    clean = np.zeros((N_STAMP, N_STAMP), dtype=np.int32)
    ones = np.ones((N_STAMP, N_STAMP))
    epochs = {f"210000{i}-1{i}": (clean.copy(), ones) for i in range(3)}
    vignet, tile_cat, psf_obj, gal_obj = _fake_inputs(epochs)
    del vignet.gal_vign_cat  # OFFSET must come from the supplied gal_obj.
    for i, name in enumerate(epochs):
        gal_obj[name]["OFFSET"] = np.array([0.1 * i, -0.1 * i])
    stamp = prepare_postage_stamps(
        vignet, 1, 0, tile_cat, bkg_sub=False,
        psf_obj=psf_obj, gal_obj=gal_obj,
    )
    names = {id(flag): name for name, (flag, _) in epochs.items()}
    assert len(stamp.offsets) == len(stamp.flags) == 3
    for flag, offset in zip(stamp.flags, stamp.offsets):
        npt.assert_array_equal(offset, gal_obj[names[id(flag)]]["OFFSET"])


# --- Ngmix.process: the per-tile epoch-cut tally ---------------------------

class _RecordingLogger:
    def __init__(self):
        self.messages = []

    def info(self, msg, *_args, **_kwargs):
        self.messages.append(msg)

    warning = error = info


def test_process_logs_the_epoch_cut_tally(tmp_path, monkeypatch):
    """One tile, four objects; the end-of-tile line counts each cut's drops.

    * object 1: clean, 18-column edge band, defect 3 px from the centre -> one
      epoch each for considered, masked_fraction, central_veto; survives.
    * object 2: edge band and central defect -> both epochs dropped; emptied.
    * object 3: one all-zero stamp, skipped before the cuts -> not considered,
      and not emptied by the cuts.
    * object 4: no PSF ('empty') -> never reaches the cuts.

    Failure modes: a cut's drops are not counted or land in the wrong
    counter; epochs skipped before the cuts are counted as considered; an
    object with no epoch at all is reported as emptied by the cuts; counts
    from one object overwrite another's.
    """
    centre = N_STAMP // 2
    clean = np.zeros((N_STAMP, N_STAMP), dtype=np.int32)
    ones = np.ones((N_STAMP, N_STAMP))
    wide18, near = clean.copy(), clean.copy()
    wide18[:, -18:] = 1
    near[centre, centre + 3] = 1
    objects = {
        1: {"2100001-10": (clean, ones), "2100002-11": (wide18, ones),
            "2100003-12": (near, ones)},
        2: {"2100004-13": (wide18.copy(), ones),
            "2100005-14": (near.copy(), ones)},
        3: {"2100006-15": (clean.copy(), ones)},
    }
    stores = {}
    for obj_id, epochs in objects.items():
        vignet, _, psf_obj, gal_obj = _fake_inputs(epochs)
        if obj_id == 3:
            gal_obj["2100006-15"]["VIGNET"] = np.zeros((N_STAMP, N_STAMP))
        stores[obj_id] = (vignet, psf_obj, gal_obj)
    vignet = SimpleNamespace(
        bkg_vign_cat=None,
        bkg_rms_vign_cat=None,
        psf_vign_cat={
            "4": "empty", **{str(i): s[1] for i, s in stores.items()}
        },
        gal_vign_cat={
            "4": "empty", **{str(i): s[2] for i, s in stores.items()}
        },
        flag_vign_cat={
            str(i): s[0].flag_vign_cat["1"] for i, s in stores.items()
        },
        weight_vign_cat={
            str(i): s[0].weight_vign_cat["1"] for i, s in stores.items()
        },
        f_wcs_file={
            k: v for s in stores.values() for k, v in s[0].f_wcs_file.items()
        },
        close=lambda: None,
    )
    tile_cat = SimpleNamespace(
        obj_id=np.array([1, 2, 3, 4]), ra=np.full(4, RA), dec=np.full(4, DEC),
        vign=None, seg=None, flux=None,
    )

    paths = [tmp_path / f"{name}.sqlite" for name in
             ("gal", "psf", "weight", "flag", "headers")]
    for path in paths:
        SqliteDict(str(path)).close()
    log = _RecordingLogger()
    ngmix = Ngmix(
        ["tile_cat.fits"] + [str(p) for p in paths[:4]],
        str(tmp_path), "-001-001", 30.0, str(paths[4]), log,
        bkg_sub=False,
    )
    ngmix._vignet_cat.close()
    ngmix._vignet_cat = vignet
    monkeypatch.setattr(ngmix_module, "Tile_cat", lambda *a, **k: tile_cat)

    def no_fit(*_args, **_kwargs):
        raise RuntimeError("metacal is not under test")

    monkeypatch.setattr(ngmix_module, "do_ngmix_metacal", no_fit)
    for method in ("compile_results", "save_results", "log_mean_ellipticity"):
        monkeypatch.setattr(Ngmix, method, lambda *_a, **_k: None)

    ngmix.process()

    lines = [m for m in log.messages if m.startswith("epoch cuts:")]
    assert len(lines) == 1, log.messages
    tally = dict(
        (k, int(v)) for k, v in re.findall(r"(\w+)=(\d+)", lines[0])
    )
    assert tally == dict(
        considered=5, masked_fraction=2, central_veto=2, objects_emptied=1
    ), lines[0]


# --- make_ngmix_observation: the HSM centroid reads the filled image -------

def test_hsm_centroid_ignores_raw_defect_values():
    """With ``centroid_source="hsm"``, the Jacobian origin lands on the
    object even when a flagged column a few pixels away holds a raw value
    ten times the object's peak.

    Failure mode: HSM measures the raw stamp, so the defect drags the centroid
    off the object (or HSM fails and falls back to the stamp centre).
    """
    import galsim

    n = 51
    centre = (n - 1) / 2
    d_row, d_col = 1.3, -0.7
    rows, cols = np.mgrid[:n, :n]
    gal = np.exp(
        -((rows - centre - d_row) ** 2 + (cols - centre - d_col) ** 2)
        / (2 * 2.0 ** 2)
    )
    flag = np.zeros((n, n), dtype=np.int32)
    flag[:, n // 2 + 6] = 1
    gal[flag != 0] = 10.0
    psf = np.exp(-((rows - centre) ** 2 + (cols - centre) ** 2) / 2.0)
    obs = make_ngmix_observation(
        gal, np.ones((n, n)), flag, psf / psf.sum(),
        galsim.PixelScale(0.1857).jacobian(), np.random.RandomState(0),
        bkg_rms=np.full((n, n), 1e-3), centroid_source="hsm",
    )
    row, col = obs.jacobian.get_cen()
    npt.assert_allclose(
        [row - centre, col - centre], [d_row, d_col], atol=0.05
    )


# --- Interpolated defects ---------------------------------------------------

def _hot_stamp():
    """A bright object with hot defects: a column and a finite 3-px bleed
    (interpolated), an edge band and a 5x5 blob (noise-filled)."""
    centre = N_STAMP // 2
    rows, cols = np.mgrid[:N_STAMP, :N_STAMP]
    gal = 1e3 * np.exp(
        -((rows - centre) ** 2 + (cols - centre) ** 2) / (2 * 3.0 ** 2)
    )
    flag = np.zeros((N_STAMP, N_STAMP), dtype=np.int32)
    flag[:, centre + 8] = 1
    flag[centre - 5:centre + 6, centre - 11:centre - 8] = 1
    flag[:, :4] = 1
    flag[centre + 9:centre + 14, centre + 12:centre + 17] = 1
    gal[flag != 0] = 5e4
    return gal, np.ones((N_STAMP, N_STAMP)), flag


def _expected_masks(target, defect, defect_weighting):
    """The fill and the zero-weight set each DEFECT_WEIGHTING promises,
    written out independently of ``defect_weighting_masks``."""
    noisefilled = defect & ~target
    if defect_weighting == "des_y6":
        return target | (np.rot90(target) & ~defect), noisefilled
    zero = {
        "fourfold_zero": defect | fourfold(target),
        "hole": defect,
        "full": noisefilled,
    }[defect_weighting]
    return target, zero


@pytest.mark.parametrize("defect_weighting", DEFECT_WEIGHTINGS)
@pytest.mark.parametrize("blend_handling", ["noisefill", "uberseg"])
def test_interpolated_fill_and_its_weights(blend_handling, defect_weighting):
    """Short defect runs take the interpolant of the clean image; other
    defects take noise. Under fourfold_zero the weight is zero on the
    defects and on the quarter-turn orbit of the interpolated pixels, whose
    light stays; under hole on the defects; under full and des_y6 only on
    the noise-filled ones. des_y6 also interpolates the clean quarter-turn
    (k=1) copies of the interpolated pixels, at full weight. The uberseg
    neighbour side stays at weight 0 under every option. The metacal noise
    image is interpolated by the same operator in fixnoise's quarter-turned
    frame.

    Failure modes: the orbit is not zero-weighted under fourfold_zero (a
    one-sided hole in the likelihood biases c), or its light is replaced;
    the copies are not interpolated under des_y6, or a noise-filled pixel
    is; an interpolated pixel on the neighbour side gets weight; wide
    defects or edge bands are extrapolated; raw defect values leak; the
    noise image is interpolated in the detector frame
    (noise-image-in-the-fixnoise-frame).
    """
    gal, weight, flag = _hot_stamp()
    defect = flag != 0
    target = interpolable_defects(defect)
    assert target.any() and (defect & ~target).any()
    kwargs = (
        dict(seg=_uberseg_seg(N_STAMP), object_number=1, dilate_neighbour=1)
        if blend_handling == "uberseg"
        else {}
    )
    gal_out, w_out, noise_out = prepare_ngmix_weights(
        gal, weight, flag, np.random.RandomState(4),
        bkg_rms=np.ones((N_STAMP, N_STAMP)), blend_handling=blend_handling,
        defect_weighting=defect_weighting, **kwargs,
    )

    neighbour = (
        uberseg_mask(kwargs["seg"], 1, dilate_neighbour=1)
        if kwargs else np.zeros_like(defect)
    )
    fill, zero = _expected_masks(target, defect, defect_weighting)
    if defect_weighting == "des_y6":
        assert (fill & ~defect).any()
    npt.assert_array_equal(w_out == 0.0, zero | neighbour)
    npt.assert_array_equal(w_out[~(zero | neighbour)], 1.0)
    kept = ~defect & ~fill
    npt.assert_array_equal(gal_out[kept], gal[kept])
    expected = interpolate_defects(gal[None], defect | fill, fill)[0]
    npt.assert_allclose(gal_out[fill], expected[fill], rtol=1e-5)
    assert np.all(np.abs(gal_out[defect & ~target]) < 10.0)
    turned = np.rot90(noise_out)
    refilled = interpolate_defects(turned[None], defect | fill, fill)[0]
    npt.assert_allclose(turned[fill], refilled[fill], atol=1e-5)


def test_fixnoise_turns_the_noise_image_k1_then_k3():
    """ngmix's fixnoise adds np.rot90(sheared(np.rot90(noise, 1)), 3) to
    each sheared image: with an asymmetric noise image, the fixnoise output
    minus the plain metacal output is exactly that.

    Failure mode: an ngmix release changes the turn (k=3 first, or none),
    and the noise image prepare_ngmix_weights interpolates in the k=1 frame
    no longer lines up with the science image's interpolated pixels
    (noise-image-in-the-fixnoise-frame).
    """
    from ngmix import DiagonalJacobian, Observation
    from ngmix.metacal import get_all_metacal

    n = 33
    rows, cols = np.indices((n, n)) - (n - 1) / 2
    psf_image = np.exp(-(rows ** 2 + cols ** 2) / (2 * 2.0 ** 2))
    psf_image /= psf_image.sum()
    galaxy = 100 * np.exp(-(rows ** 2 + cols ** 2) / (2 * 3.0 ** 2))
    jacobian = DiagonalJacobian(row=(n - 1) / 2, col=(n - 1) / 2, scale=0.2)
    noise = np.random.RandomState(3).normal(size=(n, n))
    noise[:, 20] += 5.0

    def metacal(image, noise_image=None):
        obs = Observation(
            image, weight=np.ones((n, n)), jacobian=jacobian,
            psf=Observation(psf_image, jacobian=jacobian), noise=noise_image,
        )
        return get_all_metacal(
            obs, psf="gauss", types=["noshear", "1p"],
            fixnoise=noise_image is not None, use_noise_image=True,
            rng=np.random.RandomState(0),
        )

    fixed = metacal(galaxy, noise)
    plain = metacal(galaxy)
    turned = metacal(np.rot90(noise, 1).copy())
    for t in ("noshear", "1p"):
        added = fixed[t].image - plain[t].image
        npt.assert_allclose(added, np.rot90(turned[t].image, 3), atol=1e-9)
        assert not np.allclose(added, np.rot90(turned[t].image, 1))


# --- The tile VIGNET's -1e30 neighbour markers are not defects -------------
#
# The tile VIGNET carries -1e30 on the footprints of other detections. Every
# epoch shares that tile stamp, so a marker counted as a defect would drop
# every epoch of an object with a neighbour inside the veto radius. The
# markers are their own per-epoch mask (``stamp.neighbours``): noisefill
# zero-weights and noise-fills them, uberseg ignores them, and the epoch cuts
# never read them.

_MARKER = -1.0e30
_CENTRE = N_STAMP // 2
# Flipped (ccd < 18) and unflipped (ccd >= 18) MegaCam CCDs.
_MARKER_EPOCH_NAMES = ["2100001-10", "2100002-20", "2100003-11"]


def _tile_with_neighbour(columns_from=_CENTRE + 3, rows=(_CENTRE - 2,
                                                         _CENTRE + 3)):
    """Tile VIGNET with a -1e30 neighbour footprint whose nearest pixel is
    3 px from the stamp centre. The footprint is off-centre, so the MegaCam
    flip moves it."""
    tile = np.random.default_rng(5).normal(0.0, 1.0, (N_STAMP, N_STAMP))
    tile[rows[0]:rows[1], columns_from:columns_from + 6] = _MARKER
    return tile


def _marker_stamp(tile, epochs=None, **kwargs):
    """Run prepare_postage_stamps on defect-free epochs under ``tile``."""
    clean = np.zeros((N_STAMP, N_STAMP), dtype=np.int32)
    ones = np.ones((N_STAMP, N_STAMP))
    if epochs is None:
        epochs = {name: (clean.copy(), ones) for name in _MARKER_EPOCH_NAMES}
    vignet, tile_cat, psf_obj, gal_obj = _fake_inputs(epochs)
    tile_cat.vign = tile[np.newaxis]
    stamp = prepare_postage_stamps(
        vignet, 1, 0, tile_cat, bkg_sub=False,
        psf_obj=psf_obj, gal_obj=gal_obj, **kwargs,
    )
    return stamp, epochs, gal_obj


def _expected_neighbours(tile, name):
    return Ngmix.MegaCamFlip(tile, int(name.split("-")[1])) == _MARKER


def test_neighbour_markers_near_the_centre_keep_every_epoch():
    """A neighbour footprint 3 px from the centre, and one covering 41% of
    the stamp, drop no epoch: the masked-fraction cut and the central veto
    count no marker.

    Failure mode: the markers are written into the flag stamp and counted as
    defects, so every epoch (they all share the tile VIGNET) is dropped by
    the veto or the fraction cut and the object loses its shape
    (neighbour-markers-are-not-defects).
    """
    small = _tile_with_neighbour()
    # Large, but short of the stamp border: no row or column is entirely
    # marked, so it is a neighbour, not off-tile.
    large = _tile_with_neighbour()
    large[1:-1, _CENTRE + 3:-1] = _MARKER
    assert (large == _MARKER).mean() > EPOCH_MASKED_FRACTION_CUT
    for tile in (small, large):
        stamp, _, _ = _marker_stamp(tile)
        assert len(stamp.gals) == len(_MARKER_EPOCH_NAMES)
        assert stamp.epoch_cuts["considered"] == len(_MARKER_EPOCH_NAMES)
        assert stamp.epoch_cuts["masked_fraction"] == 0
        assert stamp.epoch_cuts["central_veto"] == 0


def test_neighbour_markers_are_their_own_per_epoch_mask():
    """``stamp.neighbours`` holds the MegaCam-flipped marker mask of each
    surviving epoch, and the flag stamps stay the exposure's own.

    Failure mode: the markers are merged into the flags, or the neighbour
    mask is not flipped with its epoch and lands on the wrong pixels.
    """
    tile = _tile_with_neighbour()
    stamp, epochs, _ = _marker_stamp(tile)
    names = {id(v[0]): name for name, v in epochs.items()}
    assert len(stamp.neighbours) == len(stamp.flags)
    for flag, neighbour in zip(stamp.flags, stamp.neighbours):
        name = names[id(flag)]
        npt.assert_array_equal(flag, 0)
        npt.assert_array_equal(neighbour, _expected_neighbours(tile, name))
    assert not np.array_equal(stamp.neighbours[0], stamp.neighbours[1])


def test_a_flagged_column_near_the_centre_is_still_vetoed():
    """With a neighbour footprint present, an epoch with a genuinely flagged
    column 3 px from the centre is still dropped, and only that epoch.

    Failure mode: handling the markers apart also exempts real defects from
    the central veto.
    """
    clean = np.zeros((N_STAMP, N_STAMP), dtype=np.int32)
    ones = np.ones((N_STAMP, N_STAMP))
    column = clean.copy()
    column[:, _CENTRE - 3] = 1
    epochs = {
        "2100001-10": (clean, ones),
        "2100002-20": (column, ones),
        "2100003-11": (clean.copy(), ones),
    }
    stamp, _, _ = _marker_stamp(_tile_with_neighbour(), epochs)
    names = {id(v[0]): name for name, v in epochs.items()}
    assert sorted(names[id(f)] for f in stamp.flags) == [
        "2100001-10", "2100003-11"
    ]
    assert stamp.epoch_cuts["central_veto"] == 1
    assert stamp.epoch_cuts["masked_fraction"] == 0


def _stamp_epoch_weights(stamp, i, seed, **kwargs):
    return prepare_ngmix_weights(
        1.0e3 + stamp.gals[i], stamp.weights[i], stamp.flags[i],
        np.random.RandomState(seed), bkg_rms=stamp.bkg_rms[i],
        neighbour=stamp.neighbours[i], **kwargs,
    )


def test_noisefill_fills_exactly_the_marked_pixels():
    """Under noisefill, a defect-free epoch has zero weight and noise
    exactly on the marked neighbour pixels; every other pixel keeps its raw
    value and its weight.

    Failure mode: the neighbour markers are dropped with the flags, so
    noisefill no longer removes neighbour light (noisefill-fills-markers).
    """
    tile = _tile_with_neighbour()
    stamp, epochs, _ = _marker_stamp(tile)
    assert len(stamp.gals) == len(_MARKER_EPOCH_NAMES)
    for i in range(len(stamp.gals)):
        gal = 1.0e3 + stamp.gals[i]
        gal_out, w_out, _ = _stamp_epoch_weights(
            stamp, i, seed=i, blend_handling="noisefill",
        )
        neighbour = stamp.neighbours[i]
        assert neighbour.any()
        npt.assert_array_equal(gal_out != gal, neighbour)
        npt.assert_array_equal(w_out == 0.0, neighbour)
        assert np.all(np.abs(gal_out[neighbour]) < 10.0)


def test_uberseg_leaves_the_marked_pixels_raw():
    """Under uberseg, the markers mask nothing: with a seg map holding only
    the central object, a defect-free epoch keeps every pixel raw and
    weighted, marked or not.

    Failure mode: the markers reach the defect set or the fill, so uberseg
    noise-fills the neighbour's light instead of leaving it to the seg-based
    weight (uberseg-ignores-markers).
    """
    tile = _tile_with_neighbour()
    stamp, _, _ = _marker_stamp(tile)
    seg = np.zeros((N_STAMP, N_STAMP), dtype=np.int32)
    seg[_CENTRE - 1:_CENTRE + 2, _CENTRE - 1:_CENTRE + 2] = 1
    assert len(stamp.gals) == len(_MARKER_EPOCH_NAMES)
    for i in range(len(stamp.gals)):
        gal = 1.0e3 + stamp.gals[i]
        gal_out, w_out, _ = _stamp_epoch_weights(
            stamp, i, seed=i, blend_handling="uberseg", seg=seg,
            object_number=1,
        )
        npt.assert_array_equal(gal_out, gal)
        assert np.all(w_out > 0.0)


def test_do_ngmix_metacal_threads_each_epochs_neighbour_mask(monkeypatch):
    """Each epoch's neighbour mask reaches make_ngmix_observation.

    Failure mode: the mask is built but never used, so noisefill silently
    stops filling neighbours.
    """
    tile = _tile_with_neighbour()
    stamp, _, _ = _marker_stamp(tile)
    assert len(stamp.gals) == len(_MARKER_EPOCH_NAMES)
    seen = []

    class _Stop(Exception):
        pass

    def fake_observation(*args, **kwargs):
        seen.append(kwargs["neighbour"])
        if len(seen) == len(stamp.gals):
            raise _Stop
        return None

    monkeypatch.setattr(
        ngmix_module, "make_ngmix_observation", fake_observation,
    )
    monkeypatch.setattr(ngmix_module, "ObsList", list)
    with pytest.raises(_Stop):
        ngmix_module.do_ngmix_metacal(
            stamp, None, 1.0, np.random.RandomState(0),
        )
    for got, want in zip(seen, stamp.neighbours):
        assert got is want


# --- Off-tile pixels are defects -------------------------------------------
#
# The tile VIGNET also holds -1e30 beyond the tile's edge, where the epoch
# holds the object's own light, cut off. Those pixels are the runs of
# entirely -1e30 stamp rows and columns that start at a stamp border (the
# off-image part of a rectangle clip); they join the epoch's defect set with
# zero exposure weight. The other markers are the neighbour mask.


def _off_tile_expected(tile, name):
    flipped = Ngmix.MegaCamFlip(tile, int(name.split("-")[1])) == _MARKER
    return flipped.all(axis=1)[:, None] | flipped.all(axis=0)[None, :]


def test_object_three_px_from_the_tile_edge_is_dropped():
    """With the tile edge 3 px from the object, every epoch is dropped: the
    off-tile band counts toward the epoch cuts.

    Failure mode: off-tile pixels are treated as neighbour markers, so an
    edge object is measured with a noise-filled band through its own light
    (off-tile-pixels-are-defects).
    """
    tile = np.random.default_rng(5).normal(0.0, 1.0, (N_STAMP, N_STAMP))
    tile[:, :_CENTRE - 2] = _MARKER
    stamp, _, _ = _marker_stamp(tile)
    assert len(stamp.gals) == 0
    assert stamp.epoch_cuts["considered"] == len(_MARKER_EPOCH_NAMES)
    assert (
        stamp.epoch_cuts["masked_fraction"] + stamp.epoch_cuts["central_veto"]
        == len(_MARKER_EPOCH_NAMES)
    )


def test_the_central_veto_sees_off_tile_pixels(monkeypatch):
    """An off-tile band at EPOCH_CENTRAL_DEFECT_RADIUS from the object passes
    the central veto; one 9 px away is noise-filled inside it and vetoed.

    Every off-tile band within reach of the veto also exceeds the 10%
    fraction cut (this one is 17 columns, 1/3 of the stamp), so the fraction
    cut is lifted here to isolate the veto, which is what keeps an edge
    object's epochs safe if the fraction cut is ever relaxed.

    Failure mode: the central veto does not read the off-tile set.
    """
    monkeypatch.setattr(ngmix_module, "EPOCH_MASKED_FRACTION_CUT", 1.0)
    sky = np.random.default_rng(5).normal(0.0, 1.0, (N_STAMP, N_STAMP))
    far, near = sky.copy(), sky.copy()
    rn = int(EPOCH_CENTRAL_DEFECT_RADIUS)
    far[:, :_CENTRE - rn + 1] = _MARKER
    near[:, :_CENTRE - 8] = _MARKER
    assert (near == _MARKER).mean() > EPOCH_MASKED_FRACTION_CUT
    kept, _, _ = _marker_stamp(far)
    assert len(kept.gals) == len(_MARKER_EPOCH_NAMES)
    vetoed, _, _ = _marker_stamp(near)
    assert len(vetoed.gals) == 0
    assert vetoed.epoch_cuts["central_veto"] == len(_MARKER_EPOCH_NAMES)


def test_corner_off_tile_region_and_border_neighbour_are_classified():
    """At a tile corner, exactly the L-shaped off-tile region gets zero
    weight, and a neighbour footprint touching the stamp border without
    filling a row or column stays in the neighbour mask.

    Failure modes: off-tile pixels are classified by something other than
    whole marked rows and columns (the L is missed or a border-touching
    neighbour is swallowed); the classification ignores the MegaCam flip.
    """
    tile = np.random.default_rng(5).normal(0.0, 1.0, (N_STAMP, N_STAMP))
    tile[:2, :] = _MARKER
    tile[:, -2:] = _MARKER
    tile[40:, :6] = _MARKER  # neighbour on the bottom-left border
    stamp, epochs, _ = _marker_stamp(tile)
    names = {id(v[0]): name for name, v in epochs.items()}
    assert len(stamp.gals) == len(_MARKER_EPOCH_NAMES)
    for flag, weight, neighbour in zip(
        stamp.flags, stamp.weights, stamp.neighbours
    ):
        name = names[id(flag)]
        off_tile = _off_tile_expected(tile, name)
        assert off_tile.sum() == 2 * N_STAMP * 2 - 4
        npt.assert_array_equal(weight == 0, off_tile)
        npt.assert_array_equal(flag, 0)
        npt.assert_array_equal(
            neighbour, _expected_neighbours(tile, name) & ~off_tile
        )
        assert neighbour.sum() == 11 * 6


@pytest.mark.parametrize("blend_handling", ["noisefill", "uberseg"])
def test_off_tile_pixels_are_zero_weighted_and_filled(blend_handling):
    """Off-tile pixels are zero-weighted and noise-filled under either
    BLEND_HANDLING; under uberseg the neighbour markers stay raw.

    Failure mode: under uberseg the off-tile band keeps its weight, or the
    fill differs between blend handlings.
    """
    tile = np.random.default_rng(5).normal(0.0, 1.0, (N_STAMP, N_STAMP))
    tile[:5, :] = _MARKER
    tile[20:24, 35:40] = _MARKER
    stamp, _, _ = _marker_stamp(tile)
    seg = np.zeros((N_STAMP, N_STAMP), dtype=np.int32)
    seg[_CENTRE - 1:_CENTRE + 2, _CENTRE - 1:_CENTRE + 2] = 1
    kwargs = (
        dict(seg=seg, object_number=1) if blend_handling == "uberseg" else {}
    )
    for i in range(len(stamp.gals)):
        gal = 1.0e3 + stamp.gals[i]
        gal_out, w_out, _ = _stamp_epoch_weights(
            stamp, i, seed=i, blend_handling=blend_handling, **kwargs,
        )
        off_tile = stamp.weights[i] == 0
        assert off_tile.sum() == 5 * N_STAMP
        removed = off_tile | (
            stamp.neighbours[i] if blend_handling == "noisefill" else False
        )
        npt.assert_array_equal(gal_out != gal, removed)
        npt.assert_array_equal(w_out == 0.0, removed)


# --- Tile VIGNETs cut from a segmentation map ------------------------------
#
# SExtractor's tile VIGNET holds -1e30 off the image and on other objects'
# segmentation footprints. Cutting stamps the same way from a segmentation
# map gives the off-image pixels (off-tile defects) and the footprint pixels
# (the neighbour mask) exactly.

DR6_PATCH = Path(__file__).parent / "data" / "dr6_202.301_seg_patch.fits"


def _segmap_stamps(seg, x, y):
    """SExtractor-style tile VIGNETs on a unit image, and each stamp's
    off-image mask. Each object owns the footprint under its centre pixel."""
    half = N_STAMP // 2
    ny, nx = seg.shape
    vignets, off_image = [], []
    for xi, yi in zip(x, y):
        r0, c0 = int(np.rint(yi)) - 1, int(np.rint(xi)) - 1
        rows = r0 - half + np.arange(N_STAMP)
        cols = c0 - half + np.arange(N_STAMP)
        off = (
            ((rows < 0) | (rows >= ny))[:, None]
            | ((cols < 0) | (cols >= nx))[None, :]
        )
        labels = seg[np.clip(rows, 0, ny - 1)[:, None],
                     np.clip(cols, 0, nx - 1)[None, :]]
        own = seg[r0, c0]
        other = (labels != 0) & (labels != own)
        vignets.append(
            np.where(off | other, _MARKER, 1.0).astype(np.float32)
        )
        off_image.append(off)
    return np.array(vignets), off_image


def test_dr6_segmap_stamps_split_into_off_image_and_neighbours():
    """On the real 202.301 segmentation patch, every stamp splits into its
    off-image pixels (off-tile) and its other -1e30 pixels (neighbours).

    Failure modes: a float32 -1e30 is not recognised as a marker;
    off-image pixels of an edge stamp land in the neighbour mask;
    neighbour-footprint pixels become off-tile defects.
    """
    with fits.open(DR6_PATCH) as hdul:
        seg = hdul["SEG"].data
        objects = hdul["OBJECTS"].data
    x, y = np.array(objects["X_IMAGE"]), np.array(objects["Y_IMAGE"])
    vignets, off_image = _segmap_stamps(seg, x, y)
    assert vignets.dtype == np.float32
    n_edge = 0
    for vign, off in zip(vignets, off_image):
        neighbour, off_tile = split_tile_markers(vign, vign.shape)
        npt.assert_array_equal(off_tile, off)
        npt.assert_array_equal(neighbour, (vign == _MARKER) & ~off)
        n_edge += off.any()
    assert n_edge >= 5
    assert sum(
        split_tile_markers(v, v.shape)[0].sum() for v in vignets
    ) > 1000


def test_a_neighbour_completing_rows_beside_the_tile_edge_stays_a_neighbour():
    """An object 20 px from the tile's left edge, with a wide neighbour
    footprint that runs from the tile edge across the stamp: in the rows of
    that footprint every stamp pixel is -1e30, off the image or on the
    neighbour. Only the off-image columns are off-tile; the footprint is the
    neighbour mask, through prepare_postage_stamps.

    Failure mode: every entirely -1e30 row counts as off-tile, so the
    neighbour's rows become defects that the epoch cuts count and the defect
    fill interpolates (off-tile-is-marked-border-rows-and-columns).
    """
    seg = np.zeros((80, 80), np.int32)
    seg[20:24, 0:50] = 5
    seg[27:32, 19:24] = 1
    x, y = np.array([21.0, 30.0]), np.array([30.0, 22.0])
    vignets, off_image = _segmap_stamps(seg, x, y)
    tile, off = vignets[0], off_image[0]
    footprint = (tile == _MARKER) & ~off
    assert footprint.sum() == 4 * (N_STAMP - 5)
    assert ((tile == _MARKER).all(axis=1) & ~off.all(axis=1)).sum() == 4

    stamp, epochs, _ = _marker_stamp(tile)
    names = {id(v[0]): name for name, v in epochs.items()}
    assert len(stamp.flags) == len(_MARKER_EPOCH_NAMES)
    for flag, weight, neighbour in zip(
        stamp.flags, stamp.weights, stamp.neighbours
    ):
        ccd = int(names[id(flag)].split("-")[1])
        npt.assert_array_equal(weight == 0, Ngmix.MegaCamFlip(off, ccd))
        npt.assert_array_equal(neighbour, Ngmix.MegaCamFlip(footprint, ccd))


# --- Interpolation beside a removed neighbour -------------------------------

def _defect_beside_neighbour():
    """A single flagged pixel 6 px right of the centre, with a bright marked
    neighbour footprint starting on the next column."""
    n = 31
    c = n // 2
    gal = np.random.default_rng(3).normal(0.0, 1.0, (n, n))
    flag = np.zeros((n, n), dtype=np.int32)
    flag[c, c + 6] = 1
    neighbour = np.zeros((n, n), dtype=bool)
    neighbour[c - 2:c + 3, c + 7:c + 10] = True
    gal[neighbour] += 1.0e4
    seg = np.zeros((n, n), dtype=np.int32)
    seg[c - 1:c + 2, c - 1:c + 2] = 1
    return gal, flag, neighbour, seg, (c, c + 6)


def test_noisefill_interpolation_does_not_read_removed_neighbour_light():
    """Under noisefill, the interpolant of a defect beside a marked neighbour
    is built from the pixels the image keeps: the neighbour's light, which
    noisefill removes, is not in its support. Under uberseg the neighbour
    light is raw in the image and supports the interpolant.

    Failure mode: the support includes the removed neighbour pixels, so the
    defect is filled with light the image no longer contains
    (interpolation-support-is-the-kept-image).
    """
    gal, flag, neighbour, seg, pix = _defect_beside_neighbour()
    weight = np.ones_like(gal)
    target = np.zeros_like(neighbour)
    target[pix] = True

    out, _, _ = prepare_ngmix_weights(
        gal, weight, flag, np.random.RandomState(0),
        blend_handling="noisefill", neighbour=neighbour,
    )
    expected = interpolate_defects(gal[None], (flag != 0) | neighbour, target)
    assert out[pix] == pytest.approx(expected[0][pix])
    assert abs(out[pix]) < 100.0

    out, _, _ = prepare_ngmix_weights(
        gal, weight, flag, np.random.RandomState(0),
        blend_handling="uberseg", seg=seg, object_number=1,
        neighbour=neighbour,
    )
    expected = interpolate_defects(gal[None], flag != 0, target)
    assert out[pix] == pytest.approx(expected[0][pix])
    assert out[pix] > 1000.0


def test_a_column_beside_a_noisefill_neighbour_is_vetoed_as_noise_filled():
    """A column 8 px from the centre is interpolated, and kept, unless a
    marked neighbour lies against it. Under noisefill the neighbour's light
    is removed, so the column's rows that end on the neighbour cannot be
    interpolated; the fill noise-fills them, and the veto drops the epoch at
    the noise-fill radius. Under uberseg the neighbour's light stays and
    supports the interpolant, so the epoch is kept.

    Failure mode: the veto decides which pixels are interpolated from the
    defect mask alone, so it keeps the epoch at the 7-px radius while the
    fill noise-fills pixels 8 px from the object, inside the radius that
    noise fill needs.
    """
    column = np.zeros((N_STAMP, N_STAMP), dtype=np.int32)
    column[:, _CENTRE + 8] = 1
    # CCD 20 is not flipped, so tile and epoch share their orientation.
    epochs = {"2100001-20": (column, np.ones((N_STAMP, N_STAMP)))}
    sky = np.random.default_rng(5).normal(0.0, 1.0, (N_STAMP, N_STAMP))
    tile = sky.copy()
    tile[_CENTRE - 3:_CENTRE + 4, _CENTRE + 9:_CENTRE + 14] = _MARKER
    for tile_vign, blend_handling, kept in (
        (sky, "noisefill", 1),
        (tile, "noisefill", 0),
        (tile, "uberseg", 1),
    ):
        stamp, _, _ = _marker_stamp(
            tile_vign, epochs, blend_handling=blend_handling,
        )
        assert len(stamp.gals) == kept
        assert stamp.epoch_cuts["central_veto"] == 1 - kept

    # The fill agrees: under noisefill only the rows clear of the neighbour
    # are interpolated (their quarter turns lose weight too); the rows
    # against it are noise-filled.
    neighbour = tile == _MARKER
    defect = column != 0
    interpolated = defect & ~(neighbour[:, _CENTRE + 9][:, None])
    _, weights, _ = prepare_ngmix_weights(
        sky, np.ones_like(sky), column, np.random.RandomState(0),
        bkg_rms=np.ones_like(sky), neighbour=neighbour,
    )
    npt.assert_array_equal(
        weights == 0.0, defect | neighbour | fourfold(interpolated)
    )
