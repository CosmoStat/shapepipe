"""Defect fill and the epoch cuts (ngmix module).

A defect is a stamp pixel with a nonzero flag, zero exposure weight or an
invalid background RMS. Before metacal, :func:`prepare_ngmix_weights` gives
every defect weight 0 and fills it, whatever ``BLEND_HANDLING`` is: with noise
at its background RMS (``DEFECT_FILL = noise``, the default) or, for short
defect runs, with an interpolant of the clean pixels around them
(``DEFECT_FILL = interpolate``). Under uberseg, pixels on the neighbour side
only lose their weight, and their image values stay raw. The per-epoch cuts
in :func:`prepare_postage_stamps` act on the same defect set: the
masked-fraction cut counts it, and the central-defect veto drops an epoch
with a defect near the stamp centre, at a radius set by the defect's fill.
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

from shapepipe.modules.ngmix_package.defect_interpolation import (
    fourfold,
    interpolable_defects,
    interpolate_defects,
)
from shapepipe.modules.ngmix_package.ngmix import (
    EPOCH_CENTRAL_DEFECT_RADIUS,
    EPOCH_INTERPOLATED_DEFECT_RADIUS,
    Ngmix,
    make_ngmix_observation,
    prepare_ngmix_weights,
    prepare_postage_stamps,
    uberseg_weight,
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
def test_filled_set_is_the_defect_set(stamp, blend_handling, seed):
    """Filled pixels are exactly the defects (flag, zero weight, bad RMS).
    They carry zero weight and look like noise. Every other pixel keeps its
    raw value, under either BLEND_HANDLING.

    Failure modes:
    * a defect source (flag, zero weight, bad RMS) is left out of the fill;
    * the filled set grows beyond the defects (for example symmetrized);
    * the filled set and the zero-weight defect set differ;
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
    assert np.all(np.abs(gal_out[filled]) < 10.0), "fill is not unit noise"
    neighbour = (
        uberseg_weight(np.ones((n, n)), kwargs["seg"], 1, dilate_neighbour=1)
        == 0.0
        if blend_handling == "uberseg"
        else np.zeros((n, n), dtype=bool)
    )
    # Zero weight on defects and neighbour side; neighbour pixels that are
    # not defects keep their raw values (filled-set equality above).
    npt.assert_array_equal(w_out == 0.0, defect | neighbour)
    npt.assert_array_equal(w_out[~(defect | neighbour)], 1.0)


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
    centre (at least 17 px for width 9)."""
    band = np.zeros((N_STAMP, N_STAMP), dtype=bool)
    band[:, :width] = True
    band[:, -width:] = True
    return band


def _epochs():
    """A clean epoch and four masked ones, all defects far from the centre.

    * band10: 10 flagged columns on one side, 19.6% (39% if symmetrized);
    * flag18: 2 x 9 flagged columns, 35.3%;
    * dead18: the same columns at zero weight, no flags;
    * rms18: the same columns with an invalid background RMS.
    """
    clean = np.zeros((N_STAMP, N_STAMP), dtype=np.int32)
    ones = np.ones((N_STAMP, N_STAMP))
    band10 = clean.copy()
    band10[:, -10:] = 1
    two = _two_sided_band(9)
    dead = ones.copy()
    dead[two] = 0.0
    rms = ones.copy()
    rms[two] = np.nan
    return {
        "2100001-10": (clean, ones),
        "2100002-11": (band10, ones),
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
    """At the default 1/3 cut, the three 35% epochs are dropped, whichever
    defect source masks them, and the 20% band survives.

    Failure modes: the cut counts flags only (keeps dead18 and rms18), omits
    one defect source, or counts a symmetrized set (drops band10); in each
    case the cut disagrees with the set that is zero-weighted and filled
    (epoch-cut-on-defect-mask).
    """
    assert _surviving(_epochs()) == ["2100001-10", "2100002-11"]


def test_epoch_cut_threshold_is_the_configured_fraction():
    """At a 10% cut (the DES Y3 and Y6 value), the 19.6% band is dropped as
    well, and only the clean epoch survives.

    Failure mode: the configured threshold is ignored.
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
    the object's light enters the fit; or the boundary is inclusive, dropping
    the column at the radius (epoch-central-defect-veto).
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


def test_process_logs_the_epoch_cut_tally(tmp_path, monkeypatch):
    """One tile, four objects; the end-of-tile line counts each cut's drops.

    * object 1: clean, 18-column edge band, defect inside the radius -> one
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
    wide18, near = clean.copy(), clean.copy()
    wide18[:, -18:] = 1
    near[centre, centre + radius - 1] = 1
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


# --- ngmix_runner: the epoch-cut options reach Ngmix ------------------------

class _OptionConfig:
    """Config stub answering from a dict; absent options take the caller's
    fallback."""

    def __init__(self, options):
        self._options = {"MAG_ZP": "30.0", "ID_OBJ_MIN": "-1",
                         "ID_OBJ_MAX": "-1", **options}

    def has_option(self, _sec, key):
        return key in self._options

    def get(self, _sec, key):
        return self._options[key]

    def getexpanded(self, _sec, key):
        return self._options[key]

    def getfloat(self, _sec, key):
        return float(self._options[key])

    def getint(self, _sec, key):
        return int(self._options[key])

    def getboolean(self, _sec, key, fallback=False):
        return fallback


def test_runner_threads_the_epoch_cut_options(tmp_path, monkeypatch):
    """EPOCH_CENTRAL_DEFECT_RADIUS, EPOCH_MASKED_FRACTION_CUT, DEFECT_FILL and
    EPOCH_INTERPOLATED_DEFECT_RADIUS reach Ngmix as configured, and default
    to the module constants (and the noise fill) when absent.

    Failure mode: the runner drops or ignores an option, so a configured A/B
    arm silently runs the default cuts.
    """
    from shapepipe.modules import ngmix_runner as runner_module

    captured = []

    class _Capture:
        def __init__(self, *_args, **kwargs):
            captured.append(kwargs)

        def process(self):
            pass

    monkeypatch.setattr(runner_module, "Ngmix", _Capture)
    inputs = [str(tmp_path / f"in{i}.sqlite") for i in range(7)]
    for path in inputs:
        SqliteDict(path).close()

    for options, radius, fraction, fill, interpolated_radius in (
        ({"EPOCH_CENTRAL_DEFECT_RADIUS": "7.5",
          "EPOCH_MASKED_FRACTION_CUT": "0.1",
          "DEFECT_FILL": "interpolate",
          "EPOCH_INTERPOLATED_DEFECT_RADIUS": "5.5"},
         7.5, 0.1, "interpolate", 5.5),
        ({}, EPOCH_CENTRAL_DEFECT_RADIUS,
         ngmix_module.EPOCH_MASKED_FRACTION_CUT, "noise",
         EPOCH_INTERPOLATED_DEFECT_RADIUS),
    ):
        runner_module.ngmix_runner(
            inputs, {"output": str(tmp_path)}, "-001-001",
            _OptionConfig(options), "NGMIX_RUNNER", _RecordingLogger(),
        )
        assert captured[-1]["epoch_central_defect_radius"] == radius
        assert captured[-1]["epoch_masked_fraction_cut"] == fraction
        assert captured[-1]["defect_fill"] == fill
        assert (captured[-1]["epoch_interpolated_defect_radius"]
                == interpolated_radius)


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


# --- DEFECT_FILL = interpolate ----------------------------------------------

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


@pytest.mark.parametrize("blend_handling", ["none", "uberseg"])
def test_interpolated_fill_and_its_weights(blend_handling):
    """Short defect runs take the interpolant of the clean image; other
    defects take noise; the weight is zero on the defects and on the
    quarter-turn orbit of the interpolated pixels, whose light stays; the
    metacal noise image is interpolated with the same operator.

    Failure modes: the orbit is not zero-weighted (a one-sided hole in the
    likelihood biases c), or its light is replaced; the fill mask is
    symmetrized; wide defects or edge bands are extrapolated; raw defect
    values leak; the noise image keeps independent noise where the science
    image is smooth.
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
        defect_fill="interpolate", **kwargs,
    )

    neighbour = (
        uberseg_weight(np.ones_like(gal), kwargs["seg"], 1, dilate_neighbour=1)
        == 0.0
        if kwargs else np.zeros_like(defect)
    )
    npt.assert_array_equal(w_out == 0.0, defect | fourfold(target) | neighbour)
    npt.assert_array_equal(gal_out[~defect], gal[~defect])
    expected = interpolate_defects(gal[None], defect, target)[0]
    npt.assert_allclose(gal_out[target], expected[target], rtol=1e-5)
    assert np.all(np.abs(gal_out[defect & ~target]) < 10.0)
    refilled = interpolate_defects(noise_out[None], defect, target)[0]
    npt.assert_allclose(noise_out[target], refilled[target], atol=1e-5)


def test_noise_is_the_default_fill():
    """Leaving DEFECT_FILL unset is bit-identical to the noise fill."""
    gal, weight, flag = _hot_stamp()
    rms = np.ones_like(gal)
    default = prepare_ngmix_weights(
        gal, weight, flag, np.random.RandomState(9), bkg_rms=rms
    )
    noise = prepare_ngmix_weights(
        gal, weight, flag, np.random.RandomState(9), bkg_rms=rms,
        defect_fill="noise",
    )
    for a, b in zip(default, noise):
        npt.assert_array_equal(a, b)


def test_unknown_defect_fill_is_rejected():
    gal, weight, flag = _hot_stamp()
    with pytest.raises(ValueError, match="DEFECT_FILL"):
        prepare_ngmix_weights(
            gal, weight, flag, np.random.RandomState(0),
            defect_fill="interp",
        )


def _fill_veto_epochs():
    """A clean epoch; a column at the interpolated-defect radius; a column
    one pixel inside it; a 5-px bleed (noise-filled) at the same radius."""
    centre = N_STAMP // 2
    radius = int(EPOCH_INTERPOLATED_DEFECT_RADIUS)
    clean = np.zeros((N_STAMP, N_STAMP), dtype=np.int32)
    ones = np.ones((N_STAMP, N_STAMP))
    at, inside, wide = clean.copy(), clean.copy(), clean.copy()
    at[:, centre + radius] = 1
    inside[:, centre + radius - 1] = 1
    wide[:, centre + radius:centre + radius + 5] = 1
    return {
        "2100001-10": (clean, ones),
        "2100002-11": (at, ones),
        "2100003-12": (inside, ones),
        "2100004-13": (wide, ones),
    }


def test_the_veto_radius_follows_the_fill():
    """Under interpolation, an interpolated defect is vetoed inside
    EPOCH_INTERPOLATED_DEFECT_RADIUS and a noise-filled one inside
    EPOCH_CENTRAL_DEFECT_RADIUS. Under the noise fill every defect keeps the
    noise radius. The masked-fraction cut counts the raw defect set in
    both modes.

    Failure modes: the interpolated radius is applied to noise-filled pixels
    (a wide hole next to the object survives), or not applied at all; the
    configured radius is ignored; the boundary is inclusive; the fraction
    cut counts the zero-weight orbit.
    """
    epochs = _fill_veto_epochs()
    assert _surviving(epochs, defect_fill="interpolate") == [
        "2100001-10", "2100002-11"
    ]
    assert _surviving(epochs) == ["2100001-10"]
    assert _surviving(
        epochs, defect_fill="interpolate",
        epoch_interpolated_defect_radius=EPOCH_INTERPOLATED_DEFECT_RADIUS + 1,
    ) == ["2100001-10"]
    assert _surviving(_epochs(), defect_fill="interpolate") == [
        "2100001-10", "2100002-11"
    ]


def test_process_threads_the_defect_fill(tmp_path, monkeypatch):
    """Ngmix hands DEFECT_FILL to both the epoch cuts and metacal: under
    interpolation the column at the interpolated-defect radius survives and
    metacal is asked to interpolate; under noise it is vetoed.

    Failure mode: the option is dropped on the way to either, so the tile
    silently runs the noise fill or the noise cuts.
    """
    epochs = {k: v for k, v in _fill_veto_epochs().items()
              if k in ("2100001-10", "2100002-11")}
    vignet, tile_cat, psf_obj, _ = _fake_inputs(epochs)
    vignet.psf_vign_cat = {"1": psf_obj}
    vignet.close = lambda: None
    tile_cat.obj_id = np.array([1])
    tile_cat.flux = None
    monkeypatch.setattr(ngmix_module, "Tile_cat", lambda *a, **k: tile_cat)
    for method in ("compile_results", "save_results", "log_mean_ellipticity"):
        monkeypatch.setattr(Ngmix, method, lambda *_a, **_k: None)
    calls = []

    def record(stamp, *_args, **kwargs):
        calls.append((len(stamp.gals), kwargs.get("defect_fill")))
        raise RuntimeError("metacal is not under test")

    monkeypatch.setattr(ngmix_module, "do_ngmix_metacal", record)
    paths = [tmp_path / f"{name}.sqlite" for name in
             ("gal", "psf", "weight", "flag", "headers")]
    for path in paths:
        SqliteDict(str(path)).close()
    for fill in ("interpolate", "noise"):
        ngmix = Ngmix(
            ["tile_cat.fits"] + [str(p) for p in paths[:4]],
            str(tmp_path), "-001-001", 30.0, 0.186, str(paths[4]),
            _RecordingLogger(), bkg_sub=False, defect_fill=fill,
        )
        ngmix._vignet_cat.close()
        ngmix._vignet_cat = vignet
        ngmix.process()
    assert calls == [(2, "interpolate"), (1, "noise")]


def test_ngmix_rejects_an_unknown_defect_fill(tmp_path):
    paths = [tmp_path / f"{name}.sqlite" for name in
             ("gal", "psf", "weight", "flag", "headers")]
    for path in paths:
        SqliteDict(str(path)).close()
    with pytest.raises(ValueError, match="DEFECT_FILL"):
        Ngmix(
            ["tile_cat.fits"] + [str(p) for p in paths[:4]],
            str(tmp_path), "-001-001", 30.0, 0.186, str(paths[4]),
            _RecordingLogger(), bkg_sub=False, defect_fill="interp",
        )
