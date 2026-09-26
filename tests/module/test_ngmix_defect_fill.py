"""Defect fill and the epoch cuts (ngmix module).

A defect is a stamp pixel with a nonzero flag, zero exposure weight or an
invalid background RMS. Before metacal, :func:`prepare_ngmix_weights` ORs the
defect mask with its 90-, 180- and 270-degree rotations. It then gives that
set weight 0 and fills it with noise, whatever ``BLEND_HANDLING`` is. Under
uberseg, pixels on the neighbour side only lose their weight, and their image
values stay raw. The per-epoch cuts in :func:`prepare_postage_stamps` act on
the same symmetrized set.
"""

import re
from types import SimpleNamespace

import numpy as np
import numpy.testing as npt
import pytest
from astropy.io import fits
from astropy.wcs import WCS
from hypothesis import given
from hypothesis import strategies as st
from sqlitedict import SqliteDict

from shapepipe.modules.ngmix_package import ngmix as ngmix_module

from shapepipe.modules.ngmix_package.ngmix import (
    EPOCH_CENTRAL_DEFECT_RADIUS,
    Ngmix,
    prepare_ngmix_weights,
    prepare_postage_stamps,
    uberseg_weight,
)


# --- prepare_ngmix_weights: the filled set is the symmetrized defect set ---

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
    # Raw values far outside the unit-RMS noise, so a filled pixel is
    # unambiguous.
    gal = 1.0e3 + np.arange(n * n, dtype=float).reshape(n, n)
    return gal, weight, flag, bkg_rms


def _uberseg_seg(n):
    """Central object on the centre pixel, neighbour footprint in a corner."""
    seg = np.zeros((n, n), dtype=np.int32)
    seg[n // 2, n // 2] = 1
    seg[:2, :2] = 2
    return seg


@given(
    stamp=defect_stamps(),
    blend_handling=st.sampled_from(["none", "uberseg"]),
    seed=st.integers(0, 2**31 - 1),
)
def test_filled_set_is_the_smallest_rot90_invariant_defect_superset(
    stamp, blend_handling, seed
):
    """The filled set contains every defect, is invariant under a 90-degree
    rotation, and holds nothing else: it is the defects' 4-fold orbit. Filled
    pixels carry zero weight and look like noise; every other pixel keeps its
    raw value, under either BLEND_HANDLING.

    Failure modes:
    * the symmetrization is dropped, reduced to one rotation, missing its
      180- or 270-degree image, or replaced by a transpose;
    * a defect source (flag, zero weight, bad RMS) is left out of the mask;
    * the filled set and the zero-weight set differ;
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

    raw = (weight == 0) | (flag != 0) | ~(np.isfinite(bkg_rms) & (bkg_rms > 0))
    orbit = raw | np.rot90(raw) | np.rot90(raw, 2) | np.rot90(raw, 3)
    if orbit.all():
        return  # fully masked: no clean pixel left to set the noise level
    filled = gal_out != gal
    assert np.all(filled[raw]), "a defect pixel was not filled"
    npt.assert_array_equal(
        filled, np.rot90(filled), "filled set is not 90-degree invariant"
    )
    assert not np.any(filled & ~orbit), "a pixel outside the orbit was filled"
    assert np.all(np.abs(gal_out[filled]) < 10.0), "fill is not unit noise"
    neighbour = (
        uberseg_weight(np.ones((n, n)), kwargs["seg"], 1, dilate_neighbour=1)
        == 0.0
        if blend_handling == "uberseg"
        else np.zeros((n, n), dtype=bool)
    )
    # Zero weight on defects and neighbour side; neighbour pixels that are
    # not defects keep their raw values (filled-set checks above).
    npt.assert_array_equal(w_out == 0.0, filled | neighbour)
    npt.assert_array_equal(w_out[~(filled | neighbour)], 1.0)


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
    gal_out, w_out, _ = prepare_ngmix_weights(
        gal, weight, flag, np.random.RandomState(0), bkg_rms=np.ones((n, n)),
        blend_handling="uberseg", seg=_uberseg_seg(n), object_number=1,
    )
    assert w_out[1, 1] == 0.0 and gal_out[1, 1] != gal[1, 1]
    # A neighbour-side pixel that is not a defect: zero weight, raw value.
    assert w_out[0, 3] == 0.0 and gal_out[0, 3] == gal[0, 3]


def test_non_square_stamp_fails_loudly():
    """Rotating by 90 degrees needs a square stamp; anything else must raise.

    Failure mode: a non-square stamp is silently mis-symmetrized, or fails
    with an unrelated broadcast error.
    """
    gal = np.ones((10, 12))
    with pytest.raises(ValueError, match="square"):
        prepare_ngmix_weights(
            gal, np.ones_like(gal), np.zeros_like(gal),
            np.random.RandomState(0),
        )


# --- prepare_postage_stamps: the fraction cut counts the symmetrized set ---

N_STAMP = 51
RA, DEC = 150.0, 2.0


def _fake_inputs(epochs):
    """Minimal vignet / tile-catalogue stand-ins for prepare_postage_stamps.

    ``epochs`` maps ``"<exp>-<ccd>"`` to ``(flag, weight)``. Returns
    ``(vignet, tile_cat, psf_obj, gal_obj)``.
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

    psf_obj = per_epoch(lambda k: np.ones((N_STAMP, N_STAMP)))
    gal_obj = per_epoch(lambda k: rng.normal(0.0, 1.0, (N_STAMP, N_STAMP)))
    vignet = SimpleNamespace(
        bkg_vign_cat=None,
        bkg_rms_vign_cat=None,
        flag_vign_cat={"1": per_epoch(lambda k: epochs[k][0])},
        weight_vign_cat={"1": per_epoch(lambda k: epochs[k][1])},
        f_wcs_file={
            k.split("-")[0]: {int(k.split("-")[1]): {"WCS": wcs, "header": header}}
            for k in epochs
        },
    )
    tile_cat = SimpleNamespace(
        vign=None, seg=None, ra=np.array([RA]), dec=np.array([DEC])
    )
    return vignet, tile_cat, psf_obj, gal_obj


def _epochs():
    """One clean epoch and three masked ones, all defects far from the centre.

    Masked fractions (raw, one rotation, then 4-fold):
    * edge3: 3 flagged columns, 5.9% -> 11.4% -> 22.1%;
    * edge5: 5 flagged columns, 9.8% -> 18.6% -> 35.4%;
    * dead5: 5 zero-weight columns, no flags, 9.8% -> 18.6% -> 35.4%.
    """
    clean = np.zeros((N_STAMP, N_STAMP), dtype=np.int32)
    ones = np.ones((N_STAMP, N_STAMP))
    edge3, edge5 = clean.copy(), clean.copy()
    edge3[:, -3:] = 1
    edge5[:, -5:] = 1
    dead5 = ones.copy()
    dead5[:, :5] = 0.0
    return {
        "2100001-10": (clean, ones),
        "2100002-11": (edge3, ones),
        "2100003-12": (edge5, ones),
        "2100004-13": (clean.copy(), dead5),
    }


def _surviving(epochs, **kwargs):
    vignet, tile_cat, psf_obj, gal_obj = _fake_inputs(epochs)
    stamp = prepare_postage_stamps(
        vignet, 1, 0, tile_cat, bkg_sub=False,
        psf_obj=psf_obj, gal_obj=gal_obj, **kwargs,
    )
    names = {id(flag): name for name, (flag, _) in epochs.items()}
    return sorted(names[id(flag)] for flag in stamp.flags)


def test_epoch_cut_counts_the_symmetrized_defect_set():
    """At the default 1/3 cut, the 5-column edge bands (35% after the 4-fold
    OR) are dropped, and the 3-column band (22%) survives.

    Failure modes: the cut counts the raw mask or a single rotation (keeps
    edge5), or counts flags only (keeps dead5). In each case an epoch
    survives whose zeroed and filled area exceeds the cut
    (epoch-cut-on-symmetrized-mask).
    """
    assert _surviving(_epochs()) == ["2100001-10", "2100002-11"]


def test_epoch_cut_threshold_is_the_configured_fraction():
    """At a 10% cut (the DES Y3/Y6 value), only the clean epoch survives. The
    3-column band sits under the cut raw (5.9%) but over it symmetrized
    (22.1%).

    Failure mode: the configured threshold is ignored, or it is applied to
    the raw mask.
    """
    assert _surviving(_epochs(), epoch_masked_fraction_cut=0.1) == [
        "2100001-10"
    ]


# --- prepare_postage_stamps: the central-defect veto -----------------------

def _veto_epochs(radius):
    """A clean epoch, one with a single flagged pixel just inside
    ``radius`` of the stamp centre, and one with a flagged column exactly
    ``radius`` away. Both masked epochs are far below the fraction cut.
    """
    centre = N_STAMP // 2
    clean = np.zeros((N_STAMP, N_STAMP), dtype=np.int32)
    ones = np.ones((N_STAMP, N_STAMP))
    near, far = clean.copy(), clean.copy()
    near[centre, centre + radius - 1] = 1
    far[:, centre + radius] = 1
    return {
        "2100001-10": (clean, ones),
        "2100002-11": (near, ones),
        "2100003-12": (far, ones),
    }


def test_central_defect_vetoes_the_epoch():
    """At the default radius, a single defect pixel inside it drops the
    epoch, and a column at the radius does not.

    Failure mode: the veto is skipped, so an epoch whose filled hole overlaps
    the object's light (a per-epoch m of -5% to -35% on the calibration sim)
    enters the fit; or the boundary is inclusive, dropping the column at the
    radius (epoch-central-defect-veto).
    """
    epochs = _veto_epochs(EPOCH_CENTRAL_DEFECT_RADIUS)
    assert _surviving(epochs) == ["2100001-10", "2100003-12"]


def test_central_defect_radius_is_the_configured_value():
    """Radius 0 disables the veto; a radius one pixel larger than the far
    column's distance drops that epoch too.

    Failure mode: the configured radius is ignored.
    """
    epochs = _veto_epochs(EPOCH_CENTRAL_DEFECT_RADIUS)
    assert _surviving(epochs, epoch_central_defect_radius=0) == [
        "2100001-10", "2100002-11", "2100003-12"
    ]
    assert _surviving(
        epochs, epoch_central_defect_radius=EPOCH_CENTRAL_DEFECT_RADIUS + 1
    ) == ["2100001-10"]


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


def test_process_logs_the_epoch_cut_tally(tmp_path, monkeypatch):
    """One tile, four objects; the end-of-tile line counts each cut's drops.

    * object 1: clean, 5-column edge band, defect inside the radius -> one
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
    radius = EPOCH_CENTRAL_DEFECT_RADIUS
    centre = N_STAMP // 2
    clean = np.zeros((N_STAMP, N_STAMP), dtype=np.int32)
    ones = np.ones((N_STAMP, N_STAMP))
    edge5, near = clean.copy(), clean.copy()
    edge5[:, -5:] = 1
    near[centre, centre + radius - 1] = 1
    objects = {
        1: {"2100001-10": (clean, ones), "2100002-11": (edge5, ones),
            "2100003-12": (near, ones)},
        2: {"2100004-13": (edge5.copy(), ones), "2100005-14": (near.copy(), ones)},
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
        psf_vign_cat={"4": "empty", **{str(i): s[1] for i, s in stores.items()}},
        gal_vign_cat={"4": "empty", **{str(i): s[2] for i, s in stores.items()}},
        flag_vign_cat={str(i): s[0].flag_vign_cat["1"] for i, s in stores.items()},
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
        str(tmp_path), "-001-001", 30.0, 0.186, str(paths[4]), log,
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
