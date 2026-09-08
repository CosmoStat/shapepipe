"""PROPERTY-BASED TESTS FOR PSF EPOCH-AVERAGING + make_cat SENTINELS.

The two PSF families this module exports (the metacal reconvolution kernel via
``average_multiepoch_psf`` and the original image PSF via
``average_original_psf``) share a single averaging core,
:func:`shapepipe.modules.ngmix_package.ngmix._average_psf_fits`. These
hypothesis properties pin the contract of that core directly — a weighted mean
over the surviving (``flags == 0``) epochs — and the companion sentinel-fill
contract on the make_cat reader, which must leave an obj_id absent from the
ngmix catalogue at its type-specific sentinel.

Strategies are kept in physically sensible ranges (positive epoch weights,
positive sizes ``T``, ``|g| < 1``) so every generated input is a valid PSF-fit
result the averaging core would actually be handed in production.
"""

import os
import tempfile
from itertools import zip_longest

import numpy as np
import numpy.testing as npt
import pytest
from astropy.io import fits
from hypothesis import given, settings
from hypothesis import strategies as st

from shapepipe.modules.make_cat_package.make_cat import SaveCatalogue
from shapepipe.modules.ngmix_package.ngmix import _average_psf_fits


# --------------------------------------------------------------------------- #
# Strategies for one valid per-epoch PSF-fit result.
# --------------------------------------------------------------------------- #

# Bounded, non-degenerate floats: positive sizes, |g| < 1, strictly positive
# weights. min_magnitude keeps weights away from 0 so wsum can never vanish on
# the surviving epochs (the core raises ZeroDivisionError there by contract).
_g_component = st.floats(min_value=-0.9, max_value=0.9, allow_nan=False)
_positive_T = st.floats(min_value=1e-3, max_value=10.0, allow_nan=False)
_positive_err = st.floats(min_value=1e-6, max_value=1.0, allow_nan=False)
_weight = st.floats(
    min_value=1e-3, max_value=1e3, allow_nan=False, allow_infinity=False
)


def _result(g1, g2, t, g1_err, g2_err, t_err, flags=0):
    """One ngmix PSF-fit result dict, as the averaging core consumes it."""
    return {
        "flags": flags,
        "g": np.array([g1, g2]),
        "g_err": np.array([g1_err, g2_err]),
        "T": t,
        "T_err": t_err,
    }


@st.composite
def _good_epoch(draw):
    """A surviving (flags == 0) epoch paired with its positive weight."""
    return (
        _result(
            draw(_g_component), draw(_g_component), draw(_positive_T),
            draw(_positive_err), draw(_positive_err), draw(_positive_err),
        ),
        draw(_weight),
    )


# --------------------------------------------------------------------------- #
# (a) weighted average lies within [min, max] of the per-epoch values.
# --------------------------------------------------------------------------- #

@settings(deadline=None, max_examples=50)
@given(st.lists(_good_epoch(), min_size=1, max_size=8))
def test_average_lies_within_per_epoch_range(epochs):
    """A positive-weight weighted mean is a convex combination of its inputs.

    For every averaged component (g1, g2, T) the result must sit within the
    [min, max] envelope of the per-epoch values. This is the defining property
    of a weighted mean with strictly positive weights; it fails immediately if
    the core ever divided by the wrong weight sum, summed unweighted, or let a
    value escape the convex hull.
    """
    out = _average_psf_fits(epochs)

    g1_vals = np.array([r["g"][0] for r, _ in epochs])
    g2_vals = np.array([r["g"][1] for r, _ in epochs])
    t_vals = np.array([r["T"] for r, _ in epochs])

    # rtol absorbs only floating-point round-off, not a real bound violation.
    npt.assert_array_less(out["g_psf"][0], g1_vals.max() + 1e-9)
    npt.assert_array_less(g1_vals.min() - 1e-9, out["g_psf"][0])
    npt.assert_array_less(out["g_psf"][1], g2_vals.max() + 1e-9)
    npt.assert_array_less(g2_vals.min() - 1e-9, out["g_psf"][1])
    npt.assert_array_less(out["T_psf"], t_vals.max() + 1e-9)
    npt.assert_array_less(t_vals.min() - 1e-9, out["T_psf"])

    assert out["n_epoch"] == len(epochs)


@settings(deadline=None, max_examples=50)
@given(st.lists(_good_epoch(), min_size=1, max_size=8))
def test_average_matches_explicit_weighted_mean(epochs):
    """The core's output equals the textbook weighted mean of the survivors.

    Stronger than the [min, max] envelope: pins the exact value, so a swap of
    weighting factor or an off-by-one in the accumulation is caught.
    """
    out = _average_psf_fits(epochs)
    w = np.array([weight for _, weight in epochs])

    g = np.array([r["g"] for r, _ in epochs])
    t = np.array([r["T"] for r, _ in epochs])
    npt.assert_allclose(out["g_psf"], (g * w[:, None]).sum(0) / w.sum())
    npt.assert_allclose(out["T_psf"], (t * w).sum() / w.sum())


# --------------------------------------------------------------------------- #
# (b) epochs with flags != 0 are excluded from the average.
# --------------------------------------------------------------------------- #

@settings(deadline=None, max_examples=50)
@given(
    st.lists(_good_epoch(), min_size=1, max_size=6),
    st.lists(
        st.tuples(_weight, st.integers(min_value=1, max_value=255)),
        min_size=1,
        max_size=6,
    ),
)
def test_flagged_epochs_are_excluded(good_epochs, flagged_specs):
    """A failed-PSF epoch (flags != 0) must not enter the average at all.

    Each flagged epoch carries poisoned NaN measurement fields and a huge,
    out-of-range T — values that would wreck the mean (NaN-poison it, or shove
    it far outside the good-epoch envelope) if they leaked in. The averaged
    result and n_epoch must match a clean average over the good epochs only,
    proving the flagged ones were dropped, not merely down-weighted.
    """
    flagged = [
        (
            _result(
                np.nan, np.nan, 1e6, np.nan, np.nan, np.nan, flags=flags
            ),
            weight,
        )
        for weight, flags in flagged_specs
    ]
    # Interleave good and flagged epochs WITHOUT duplicating any good epoch, so
    # exclusion can't be a happy accident of ordering. itertools.zip_longest
    # threads them together; the leftover tail of the longer list follows.
    mixed = [
        e
        for pair in zip_longest(good_epochs, flagged)
        for e in pair
        if e is not None
    ]

    out = _average_psf_fits(mixed)
    expected = _average_psf_fits(good_epochs)

    npt.assert_allclose(out["g_psf"], expected["g_psf"])
    npt.assert_allclose(out["T_psf"], expected["T_psf"])
    npt.assert_allclose(out["T_psf_err"], expected["T_psf_err"])
    assert out["n_epoch"] == len(good_epochs)
    assert np.isfinite(out["g_psf"]).all() and np.isfinite(out["T_psf"])


# --------------------------------------------------------------------------- #
# (c) a single surviving epoch returns that epoch's value.
# --------------------------------------------------------------------------- #

@settings(deadline=None, max_examples=50)
@given(
    _good_epoch(),
    st.lists(
        st.tuples(_weight, st.integers(min_value=1, max_value=255)),
        max_size=5,
    ),
)
def test_single_survivor_returns_its_own_value(survivor, flagged_specs):
    """When exactly one epoch survives, its values pass through untouched.

    The lone survivor's weight cancels in mean = (v*w)/w, so the result must be
    the survivor's value exactly, regardless of how many flagged epochs (with
    arbitrary weights) surround it. n_epoch must be 1.
    """
    result, weight = survivor
    flagged = [
        (_result(np.nan, np.nan, 1e6, np.nan, np.nan, np.nan, flags=f), w)
        for w, f in flagged_specs
    ]
    out = _average_psf_fits(flagged + [survivor] + flagged)

    npt.assert_allclose(out["g_psf"], result["g"])
    npt.assert_allclose(out["g_psf_err"], result["g_err"])
    npt.assert_allclose(out["T_psf"], result["T"])
    npt.assert_allclose(out["T_psf_err"], result["T_err"])
    assert out["n_epoch"] == 1


@given(
    st.lists(
        st.tuples(_weight, st.integers(min_value=1, max_value=255)),
        min_size=1, max_size=5,
    )
)
@settings(deadline=None, max_examples=25)
def test_all_epochs_failed_raises_zero_division(flagged_specs):
    """All-epochs-failed is the contract that wsum == 0 raises, not returns 0/NaN.

    An object whose PSF fit failed on every epoch reaches both PSF-family
    wrappers (reconv and orig) with zero survivors, so ``wsum`` stays 0. The
    core must raise ``ZeroDivisionError`` rather than silently divide — a
    regression returning NaN/0 would be caught here. Every input epoch is
    flagged (``flags != 0``) with an arbitrary positive weight, so no survivor
    can rescue the sum.
    """
    flagged = [
        (_result(np.nan, np.nan, 1e6, np.nan, np.nan, np.nan, flags=f), w)
        for w, f in flagged_specs
    ]
    with pytest.raises(ZeroDivisionError):
        _average_psf_fits(flagged)


# --------------------------------------------------------------------------- #
# (d) make_cat: an obj_id absent from the ngmix cat keeps its sentinel fill.
# --------------------------------------------------------------------------- #

# The sentinel value each column family is pre-filled with before the matched
# rows are overwritten (mirrors make_cat._save_ngmix_data). The property: a row
# whose obj_id never appears among the ngmix ids keeps exactly these.
_SENTINELS = {
    "NGMIX_T_NOSHEAR": 0.0,
    "NGMIX_SNR_NOSHEAR": 0.0,
    "NGMIX_FLAGS_NOSHEAR": 0.0,
    "NGMIX_T_PSF_ORIG_NOSHEAR": 0.0,
    "NGMIX_T_PSF_RECONV_NOSHEAR": 0.0,
    "NGMIX_FLUX_ERR_NOSHEAR": -1.0,
    "NGMIX_MAG_ERR_NOSHEAR": -1.0,
    "NGMIX_G1_NOSHEAR": -10.0,
    "NGMIX_G2_NOSHEAR": -10.0,
    "NGMIX_G1_PSF_ORIG_NOSHEAR": -10.0,
    "NGMIX_G1_PSF_RECONV_NOSHEAR": -10.0,
    "NGMIX_T_ERR_NOSHEAR": 1e30,
    "NGMIX_T_ERR_PSF_ORIG_NOSHEAR": 1e30,
    "NGMIX_T_ERR_PSF_RECONV_NOSHEAR": 1e30,
    "NGMIX_N_EPOCH": 0.0,
    "NGMIX_MCAL_FLAGS": 0.0,
    "NGMIX_MCAL_TYPES_FAIL": 0.0,
    "NGMIX_NEIGHBOUR_FLAG": 0.0,
}

# Per-key write format and a measured value distinct from every sentinel, so a
# matched row is unmistakably "overwritten" and an absent row unmistakably not.
_NGMIX_KEYS = [
    "id", "n_epoch_model", "mcal_types_fail", "neighbour_flag", "nfev_fit",
    "g1", "g1_err", "g2", "g2_err", "T", "T_err",
    "flux", "flux_err", "s2n", "mag", "mag_err", "flags", "mcal_flags",
    "g1_psf_orig", "g2_psf_orig", "g1_err_psf_orig", "g2_err_psf_orig",
    "T_psf_orig", "T_err_psf_orig",
    "g1_psf_reconv", "g2_psf_reconv", "g1_err_psf_reconv", "g2_err_psf_reconv",
    "T_psf_reconv", "T_err_psf_reconv",
]
_INT_KEYS = {
    "id", "n_epoch_model", "mcal_types_fail", "neighbour_flag", "nfev_fit",
    "flags", "mcal_flags"
}
_SHEAR_EXTS = ["1M", "1P", "2M", "2P", "NOSHEAR"]


class _NullLogger:
    def info(self, *_args, **_kwargs):
        pass


def _measured_row(obj_id):
    """One fit object whose every value is far from any sentinel (5 / 0.5)."""
    row = {key: (5 if key in _INT_KEYS else 0.5) for key in _NGMIX_KEYS}
    row["id"] = obj_id
    return row


def _write_ngmix_cat(path, obj_ids):
    rows = [_measured_row(oid) for oid in obj_ids]
    hdus = [fits.PrimaryHDU()]
    for ext in _SHEAR_EXTS:
        cols = [
            fits.Column(
                name=key,
                format="K" if key in _INT_KEYS else "D",
                array=np.array([row[key] for row in rows]),
            )
            for key in _NGMIX_KEYS
        ]
        hdus.append(fits.BinTableHDU.from_columns(cols, name=ext))
    fits.HDUList(hdus).writeto(path, overwrite=True)


def _run_save_ngmix(ngmix_path, obj_id):
    inst = object.__new__(SaveCatalogue)
    inst._obj_id = np.asarray(obj_id)
    inst._output_dict = {}
    inst._cat_size_target = len(inst._obj_id)
    inst._w_log = _NullLogger()
    err_msg = inst._save_ngmix_data(str(ngmix_path))
    assert err_msg is None
    return inst._output_dict


# Distinct positive integer obj_ids; split into "fit by ngmix" vs "absent".
_distinct_ids = st.lists(
    st.integers(min_value=1, max_value=10_000),
    min_size=2,
    max_size=8,
    unique=True,
)


@settings(deadline=None, max_examples=25)
@given(_distinct_ids, st.data())
def test_absent_obj_id_keeps_sentinel_fill(all_ids, data):
    """An obj_id SExtractor saw but ngmix never fit keeps every sentinel.

    The final catalogue carries all of ``all_ids``; ngmix fit only a non-empty
    proper subset. The unfit rows must retain the exact per-column sentinel
    (0 / -10 / 1e30 / -1), while the fit rows must have been overwritten to the
    measured value — so this is not vacuously satisfied by an all-sentinel cat.
    """
    fit_ids = data.draw(
        st.lists(st.sampled_from(all_ids), min_size=1, unique=True).filter(
            lambda s: 0 < len(s) < len(all_ids)
        )
    )
    absent_ids = [oid for oid in all_ids if oid not in fit_ids]

    with tempfile.TemporaryDirectory() as tmp:
        ngmix_path = os.path.join(tmp, "ngmix.fits")
        _write_ngmix_cat(ngmix_path, fit_ids)
        out = _run_save_ngmix(ngmix_path, all_ids)

    absent_idx = [all_ids.index(oid) for oid in absent_ids]
    fit_idx = [all_ids.index(oid) for oid in fit_ids]

    for col, sentinel in _SENTINELS.items():
        arr = np.asarray(out[col])
        npt.assert_allclose(
            arr[absent_idx], sentinel,
            err_msg=f"{col}: absent rows lost their sentinel {sentinel}",
        )
        # The fit rows were overwritten — measured value 5 (int cols) or 0.5,
        # both distinct from every sentinel, so the fill wasn't global.
        assert not np.any(np.isclose(arr[fit_idx], sentinel)), (
            f"{col}: a fit row still carries the sentinel {sentinel}"
        )
