"""UNIT TESTS FOR MODULE PACKAGE: MAKE_CAT.

Drives ``SaveCatalogue._save_ngmix_data`` against a synthetic ngmix catalogue
to lock in the column grammar ``ESTIMATOR_COMPONENT[_ERR]_OBJECT[_SHEAR]``
(shapepipe#749, #761): galaxy is the implicit default object and carries no
``GAL`` token (``NGMIX_G1_NOSHEAR``, never ``NGMIX_G1_GAL_NOSHEAR``), while
the PSF families keep an explicit object token —
``NGMIX_<COMPONENT>[_ERR]_<OBJECT>_<SHEAR>`` for ``PSF_ORIG``/``PSF_RECONV`` —
plus the four OBJECT/SHEAR-less metadata columns (``NGMIX[m]_MCAL_FLAGS``,
``NGMIX_N_EPOCH``, ``NGMIX_MCAL_TYPES_FAIL``, ``NGMIX_NEIGHBOUR_FLAG``). The
original image PSF
(``PSF_ORIG``) and the metacal reconvolution kernel (``PSF_RECONV``) are
independent fits of *different* PSFs, no longer the single aliased value of
the pre-#749 code.
"""

import functools
import operator

import h5py
import numpy as np
import numpy.testing as npt
import pytest
from astropy.io import fits
from ngmix.flags import LM_FUNC_NOTFINITE, NAME_MAP
from sqlitedict import SqliteDict

from shapepipe.modules.make_cat_package import make_cat
from shapepipe.modules.make_cat_package.make_cat import SaveCatalogue
from shapepipe.modules.make_cat_runner import make_cat_runner
from shapepipe.modules.ngmix_package.ngmix import Ngmix
from shapepipe.pipeline.config import CustomParser
from shapepipe.utilities import cfis


class _NullLogger:
    def info(self, *_args, **_kwargs):
        pass

    def warning(self, *_args, **_kwargs):
        pass


class _CaptureLogger(_NullLogger):
    def __init__(self):
        self.warnings = []

    def warning(self, message, *_args, **_kwargs):
        self.warnings.append(message)


# Per-object key set that compile_results emits into every shear-type
# extension (mirrors ngmix.compile_results ``names2``). Sentinel values are
# distinct per key so a mis-routed column is caught by value, and the orig vs
# reconv PSF families carry *different* sentinels so any aliasing surfaces.
NGMIX_KEYS = [
    "id",
    "n_epoch_model",
    "n_epoch_failed",
    "mcal_types_fail",
    "neighbour_flag",
    "n_epoch_interp", "min_dist_interp", "min_dist_noisefill",
    "nfev_fit",
    "g1", "g1_err", "g2", "g2_err",
    "T", "T_err",
    "flux", "flux_err", "s2n", "mag", "mag_err",
    "flags", "mcal_flags",
    "g1_psf_orig", "g2_psf_orig",
    "g1_err_psf_orig", "g2_err_psf_orig",
    "T_psf_orig", "T_err_psf_orig",
    "g1_psf_reconv", "g2_psf_reconv",
    "g1_err_psf_reconv", "g2_err_psf_reconv",
    "T_psf_reconv", "T_err_psf_reconv",
]

SHEAR_EXTS = ["1M", "1P", "2M", "2P", "NOSHEAR"]
SHEAR_EXTS_LOWER = ["1m", "1p", "2m", "2p", "noshear"]


def _ngmix_row(obj_id):
    """One object's per-key values, distinct for orig vs reconv PSF families.

    Built so every PSF column carries a sentinel that names its source: the
    orig family from ``*_orig`` integers, the reconv family from ``*_reconv``
    integers, deliberately disjoint.
    """
    return {
        "id": obj_id,
        "n_epoch_model": 3,
        "n_epoch_failed": 1,
        "mcal_types_fail": 0,
        "neighbour_flag": 1,
        "n_epoch_interp": 2, "min_dist_interp": 7.5,
        "min_dist_noisefill": 14.25,
        "nfev_fit": 7,
        "g1": 0.10, "g1_err": 0.011, "g2": -0.20, "g2_err": 0.022,
        "T": 0.30, "T_err": 0.033,
        "flux": 100.0, "flux_err": 1.0, "s2n": 55.0,
        "mag": 21.0, "mag_err": 0.05,
        "flags": 0, "mcal_flags": 0,
        # original image PSF — its own ellipticity and size
        "g1_psf_orig": 0.0041, "g2_psf_orig": -0.0031,
        "g1_err_psf_orig": 2e-5, "g2_err_psf_orig": 3e-5,
        "T_psf_orig": 0.071, "T_err_psf_orig": 8e-4,
        # metacal reconvolution kernel — round and enlarged, distinct values
        "g1_psf_reconv": 0.0012, "g2_psf_reconv": 0.0022,
        "g1_err_psf_reconv": 1e-5, "g2_err_psf_reconv": 1e-5,
        "T_psf_reconv": 0.092, "T_err_psf_reconv": 1e-3,
    }


def _write_ngmix_cat(path, obj_ids):
    """Write a synthetic ngmix FITS with the five shear-type extensions.

    Each extension carries the exact ``compile_results`` key set, with the
    same per-object values across extensions (the PSF families are
    object-level, not per-shear).
    """
    rows = [_ngmix_row(oid) for oid in obj_ids]
    hdus = [fits.PrimaryHDU()]
    for ext in SHEAR_EXTS:
        cols = [
            fits.Column(
                name=key,
                format="K" if key in ("id", "n_epoch_model", "n_epoch_failed",
                                       "mcal_types_fail",
                                       "n_epoch_interp", "nfev_fit", "flags",
                                       "mcal_flags") else "D",
                array=np.array([row[key] for row in rows]),
            )
            for key in NGMIX_KEYS
        ]
        hdus.append(fits.BinTableHDU.from_columns(cols, name=ext))
    fits.HDUList(hdus).writeto(path, overwrite=True)


def _run_save_ngmix(ngmix_path, obj_id, cat_size_target=None, w_log=None):
    """Drive ``_save_ngmix_data`` and return its populated output dict."""
    inst = object.__new__(SaveCatalogue)
    inst._obj_id = np.asarray(obj_id)
    inst._output_dict = {}
    inst._cat_size_target = (
        len(inst._obj_id) if cat_size_target is None else cat_size_target
    )
    inst._w_log = w_log or _NullLogger()

    err_msg = inst._save_ngmix_data(str(ngmix_path))
    assert err_msg is None
    return inst._output_dict


def test_save_ngmix_data_uses_new_grammar_and_no_old_names(tmp_path):
    """Produced columns follow the new grammar; no old token survives.

    The renamed grammar drops ``_PSFo``, ``_Tpsf`` and the ellipticity
    2-vector ``NGMIX_ELL_*`` entirely (shapepipe#749). A regression to any of
    those would be a silent break for the shipped param-file consumer
    (``create_final_cat.py``), so the absence is asserted directly.
    """
    ngmix_path = tmp_path / "ngmix-0.fits"
    obj_ids = [11, 22, 33]
    _write_ngmix_cat(ngmix_path, obj_ids)

    out = _run_save_ngmix(ngmix_path, obj_ids)

    # Every per-shear family is present under the new grammar.
    for shear in SHEAR_EXTS:
        for col in (
            f"NGMIX_G1_{shear}", f"NGMIX_G2_{shear}",
            f"NGMIX_G1_ERR_{shear}", f"NGMIX_G2_ERR_{shear}",
            f"NGMIX_T_{shear}", f"NGMIX_T_ERR_{shear}",
            f"NGMIX_SNR_{shear}",
            f"NGMIX_FLUX_{shear}", f"NGMIX_FLUX_ERR_{shear}",
            f"NGMIX_MAG_{shear}", f"NGMIX_MAG_ERR_{shear}",
            f"NGMIX_FLAGS_{shear}",
            f"NGMIX_G1_PSF_ORIG_{shear}", f"NGMIX_G2_PSF_ORIG_{shear}",
            f"NGMIX_T_PSF_ORIG_{shear}",
            f"NGMIX_G1_PSF_RECONV_{shear}", f"NGMIX_G2_PSF_RECONV_{shear}",
            f"NGMIX_T_PSF_RECONV_{shear}",
        ):
            assert col in out, f"missing {col}"

    # Object-level metadata columns carry no OBJECT/SHEAR token.
    for col in (
        "NGMIX_MCAL_FLAGS", "NGMIX_N_EPOCH", "NGMIX_MCAL_TYPES_FAIL",
        "NGMIX_NEIGHBOUR_FLAG", "NGMIX_N_EPOCH_FAILED", "NGMIX_N_EPOCH_INTERP",
        "NGMIX_MIN_DIST_INTERP", "NGMIX_MIN_DIST_NOISEFILL",
    ):
        assert col in out, f"missing {col}"

    # No old name survives anywhere in the produced columns.
    for col in out:
        assert "_PSFo" not in col, col
        assert "_Tpsf" not in col and "TPSF" not in col, col
        assert not col.startswith("NGMIX_ELL_"), col
        assert "NGMIXm" not in col, col  # moments branch off by default
        # Galaxy is the implicit default object — no GAL token (shapepipe#761).
        assert "_GAL_" not in col and not col.endswith("_GAL"), col


def test_save_ngmix_data_psf_families_trace_distinct_sources(tmp_path):
    """The two PSF families read from their own ngmix columns, un-aliased.

    ``NGMIX_G1_PSF_ORIG_*`` must equal the orig source and
    ``NGMIX_G1_PSF_RECONV_*`` the reconv source, and the two must differ —
    the un-aliasing that shapepipe#749 restores. Likewise the sizes:
    ``NGMIX_T_PSF_ORIG_*`` != ``NGMIX_T_PSF_RECONV_*``.
    """
    ngmix_path = tmp_path / "ngmix-1.fits"
    obj_ids = [11, 22, 33]
    _write_ngmix_cat(ngmix_path, obj_ids)
    row = _ngmix_row(11)

    out = _run_save_ngmix(ngmix_path, obj_ids)

    for shear in SHEAR_EXTS:
        npt.assert_allclose(
            out[f"NGMIX_G1_PSF_ORIG_{shear}"],
            [row["g1_psf_orig"]] * len(obj_ids),
        )
        npt.assert_allclose(
            out[f"NGMIX_G1_PSF_RECONV_{shear}"],
            [row["g1_psf_reconv"]] * len(obj_ids),
        )
        # the un-aliasing: orig and reconv ellipticities differ
        assert not np.allclose(
            out[f"NGMIX_G1_PSF_ORIG_{shear}"],
            out[f"NGMIX_G1_PSF_RECONV_{shear}"],
        )
        # sizes are independent too
        npt.assert_allclose(
            out[f"NGMIX_T_PSF_ORIG_{shear}"], [row["T_psf_orig"]] * len(obj_ids)
        )
        npt.assert_allclose(
            out[f"NGMIX_T_PSF_RECONV_{shear}"],
            [row["T_psf_reconv"]] * len(obj_ids),
        )
        assert not np.allclose(
            out[f"NGMIX_T_PSF_ORIG_{shear}"],
            out[f"NGMIX_T_PSF_RECONV_{shear}"],
        )


def test_save_ngmix_data_fills_sentinels_for_absent_objects(tmp_path):
    """An obj_id absent from the ngmix cat keeps its sentinel fill.

    make_cat pre-fills every column with a type-specific sentinel and only
    overwrites the rows whose ``NUMBER`` matches an ngmix ``id``. An object
    SExtractor saw but ngmix never fit (no matching id) must therefore keep
    the sentinels: 0 for sizes, -10 for ellipticities, 1e30 for ``T_ERR``,
    -1 for flux/mag errors. (Its flag columns are pinned by
    ``test_galaxy_cut_admits_only_measured_objects``.)
    """
    ngmix_path = tmp_path / "ngmix-2.fits"
    # ngmix fit only object 22; the final cat also carries 11 and 99.
    _write_ngmix_cat(ngmix_path, [22])
    obj_ids = [11, 22, 99]

    out = _run_save_ngmix(ngmix_path, obj_ids, cat_size_target=3)

    present = obj_ids.index(22)
    absent = [obj_ids.index(11), obj_ids.index(99)]

    # The matched object got its measured value; the others keep sentinels.
    row = _ngmix_row(22)
    g1_orig = np.asarray(out["NGMIX_G1_PSF_ORIG_NOSHEAR"])
    npt.assert_allclose(g1_orig[present], row["g1_psf_orig"])
    npt.assert_allclose(g1_orig[absent], [-10.0, -10.0])

    t_orig = np.asarray(out["NGMIX_T_PSF_ORIG_NOSHEAR"])
    npt.assert_allclose(t_orig[present], row["T_psf_orig"])
    npt.assert_allclose(t_orig[absent], [0.0, 0.0])

    t_err_gal = np.asarray(out["NGMIX_T_ERR_NOSHEAR"])
    npt.assert_allclose(t_err_gal[present], row["T_err"])
    npt.assert_allclose(t_err_gal[absent], [1e30, 1e30])

    flux_err = np.asarray(out["NGMIX_FLUX_ERR_NOSHEAR"])
    npt.assert_allclose(flux_err[present], row["flux_err"])
    npt.assert_allclose(flux_err[absent], [-1.0, -1.0])

    n_epoch = np.asarray(out["NGMIX_N_EPOCH"])
    npt.assert_allclose(n_epoch[present], row["n_epoch_model"])
    npt.assert_allclose(n_epoch[absent], [0.0, 0.0])

    for col, key, never_fit in (
        ("NGMIX_N_EPOCH_FAILED", "n_epoch_failed", 0.0),
        ("NGMIX_N_EPOCH_INTERP", "n_epoch_interp", 0.0),
        ("NGMIX_MIN_DIST_INTERP", "min_dist_interp", -1.0),
        ("NGMIX_MIN_DIST_NOISEFILL", "min_dist_noisefill", -1.0),
    ):
        values = np.asarray(out[col])
        npt.assert_allclose(values[present], row[key])
        npt.assert_allclose(values[absent], [never_fit, never_fit])


def test_low_match_fraction_warns_and_continues_with_sentinels(tmp_path):
    """A low match count warns while unmatched detections stay in the output."""
    ngmix_path = tmp_path / "ngmix-low-match.fits"
    obj_ids = list(range(1, 12))
    _write_ngmix_cat(ngmix_path, [obj_ids[0]])
    logger = _CaptureLogger()

    out = _run_save_ngmix(
        ngmix_path,
        obj_ids,
        cat_size_target=len(obj_ids),
        w_log=logger,
    )

    assert len(logger.warnings) == 1
    assert "continuing with sentinels" in logger.warnings[0]
    assert out["NGMIX_G1_NOSHEAR"][1] == -10.0


def _metacal_result(obj_id):
    """One clean object's result as ``Ngmix.process`` hands it on.

    Every metacal type carries a successful fit; the PSF families are
    object-level, copied from ``_ngmix_row``.
    """
    row = _ngmix_row(obj_id)
    fit = {
        "nfev": row["nfev_fit"],
        "g": [row["g1"], row["g2"]],
        "g_cov": np.diag([row["g1_err"] ** 2, row["g2_err"] ** 2]),
        "T": row["T"], "T_err": row["T_err"],
        "flux": row["flux"], "flux_err": row["flux_err"],
        "s2n": row["s2n"], "flags": 0,
    }
    res = {
        "obj_id": obj_id,
        "n_epoch_model": row["n_epoch_model"],
        "n_epoch_failed": row["n_epoch_failed"],
        "neighbour_flag": row["neighbour_flag"],
        "n_epoch_interp": row["n_epoch_interp"],
        "min_dist_interp": row["min_dist_interp"],
        "min_dist_noisefill": row["min_dist_noisefill"],
    }
    for key in NGMIX_KEYS:
        if key.endswith("_psf_orig") or key.endswith("_psf_reconv"):
            res[key] = row[key]
    res.update({name: dict(fit) for name in SHEAR_EXTS_LOWER})
    return res


def _serialise_then_merge(tmp_path, results, cat_ids):
    """Write ``results`` with ngmix's own write path, read via make_cat.

    ``compile_results`` + ``save_results`` produce the ngmix catalogue;
    ``_save_ngmix_data`` merges it into a final catalogue of ``cat_ids``.
    """
    ngmix_inst = object.__new__(Ngmix)
    ngmix_inst._zero_point = 30.0
    ngmix_inst._output_dir = str(tmp_path)
    ngmix_inst._file_number_string = "-0"
    ngmix_inst.save_results(ngmix_inst.compile_results(results))
    return _run_save_ngmix(ngmix_inst.get_output_path(str(tmp_path)), cat_ids)


def test_save_ngmix_data_matches_module_serialised_catalogue(tmp_path):
    """End-to-end: make_cat reads a catalogue ngmix itself serialised.

    Rather than hand-rolling the FITS, drive ngmix's own
    ``compile_results`` + ``save_results`` (the real write path), then read it
    back through ``_save_ngmix_data``. This guards the make_cat reader against
    drift in the ngmix key set: the keys must line up end to end.
    """
    obj_ids = [5, 7]
    out = _serialise_then_merge(
        tmp_path, [_metacal_result(oid) for oid in obj_ids], obj_ids
    )

    # Both PSF families survive the round trip, distinct, on every object.
    npt.assert_allclose(
        out["NGMIX_G1_PSF_ORIG_NOSHEAR"], [_ngmix_row(o)["g1_psf_orig"] for o in obj_ids]
    )
    npt.assert_allclose(
        out["NGMIX_G1_PSF_RECONV_NOSHEAR"],
        [_ngmix_row(o)["g1_psf_reconv"] for o in obj_ids],
    )
    assert not np.allclose(
        out["NGMIX_G1_PSF_ORIG_NOSHEAR"], out["NGMIX_G1_PSF_RECONV_NOSHEAR"]
    )


# --- _save_psf_data: EXP_ID_n / CCD_n alignment with HSM_*_PSF_n (#890) ---


def _psf_epoch(g1, g2, t, flag=0, m4=None):
    """One epoch's interpolated-PSF HSM shape entry (psfex_interp's SHAPES dict).

    ``m4`` is an optional ``(M4_1, M4_2, RHO4)`` triple; omitted, the entry
    mimics a producer that predates the fourth-moment columns.
    """
    shapes = {
        "HSM_G1_PSF": g1,
        "HSM_G2_PSF": g2,
        "HSM_T_PSF": t,
        "HSM_FLAG_PSF": flag,
    }
    if m4 is not None:
        shapes.update(zip(("HSM_M4_1_PSF", "HSM_M4_2_PSF", "HSM_RHO4_PSF"), m4))
    return {"SHAPES": shapes}


def _write_galaxy_psf_cat(path, per_obj):
    """Write a synthetic ``galaxy_psf`` sqlite catalogue.

    ``per_obj`` maps object id -> ``"empty"`` or an ordered mapping of
    ``"exp_name-ccd_n"`` -> epoch entry (see ``_psf_epoch``); the mapping's
    insertion order is what defines "epoch n" (shapepipe#890).
    """
    db = SqliteDict(str(path))
    for obj_id, value in per_obj.items():
        db[str(obj_id)] = value
    db.commit()
    db.close()


def _run_save_psf(
    galaxy_psf_path, obj_id, n_overlap, n_epoch_slots=None, epoch_slots=True
):
    """Drive ``_save_psf_data`` and return its populated output dict."""
    inst = object.__new__(SaveCatalogue)
    inst._obj_id = np.asarray(obj_id)
    inst._output_dict = {}
    inst._final_cat = {"N_EPOCH": np.asarray(n_overlap)}

    inst._save_psf_data(
        str(galaxy_psf_path),
        n_epoch_slots=n_epoch_slots,
        epoch_slots=epoch_slots,
    )
    return inst._output_dict


def test_save_psf_data_exp_id_ccd_align_with_hsm_psf_slots(tmp_path):
    """EXP_ID_n/CCD_n name the same epoch as HSM_*_PSF_n, slot for slot.

    Two epochs in a known, non-alphabetical order for one object; the
    ordering of ``EXP_ID_n``/``CCD_n`` must follow the same ``key``
    iteration that assigns the HSM slots -- both come from the exact same
    enumeration in ``_save_psf_data`` (shapepipe#890), so a mis-ordering
    here would mean the two column families were built from different
    loops, not the same one.
    """
    galaxy_psf_path = tmp_path / "galaxy_psf.sqlite"
    per_obj = {
        101: {
            # Insertion order is deliberately not sorted -- an
            # index-only (not key-preserving) implementation would not
            # reproduce it by accident.
            "2229900-13": _psf_epoch(0.03, 0.04, 0.6),
            "2113864-7": _psf_epoch(0.01, 0.02, 0.5),
        },
    }
    _write_galaxy_psf_cat(galaxy_psf_path, per_obj)

    out = _run_save_psf(galaxy_psf_path, [101], n_overlap=[2])

    npt.assert_allclose(out["HSM_G1_PSF_1"], [0.03])
    npt.assert_allclose(out["HSM_G1_PSF_2"], [0.01])
    assert out["EXP_ID_1"][0] == 2229900
    assert out["CCD_1"][0] == 13
    assert out["EXP_ID_2"][0] == 2113864
    assert out["CCD_2"][0] == 7


def test_save_psf_data_identity_survives_failed_hsm_fit(tmp_path):
    """A flagged epoch keeps its EXP_ID_n/CCD_n though HSM_*_PSF_n stays sentinel.

    ``_save_psf_data`` skips writing ``HSM_*_PSF_n`` when
    ``HSM_FLAG_PSF != 0`` for that epoch (the interpolated PSF's shape fit
    failed), leaving the column at its sentinel. The epoch identity is a
    fact about which exposure/CCD occupies that slot, independent of
    whether the shape fit converged there, so EXP_ID_n/CCD_n are written
    regardless.
    """
    galaxy_psf_path = tmp_path / "galaxy_psf.sqlite"
    per_obj = {
        202: {
            "2113864-9": _psf_epoch(0.05, 0.06, 0.7),
            "2358123-21": _psf_epoch(-10.0, -10.0, 0.0, flag=5),
        },
    }
    _write_galaxy_psf_cat(galaxy_psf_path, per_obj)

    out = _run_save_psf(galaxy_psf_path, [202], n_overlap=[2])

    npt.assert_allclose(out["HSM_G1_PSF_1"], [0.05])
    # Sentinel: the epoch-2 HSM fit failed, so the pre-fill value stands.
    npt.assert_allclose(out["HSM_G1_PSF_2"], [-10.0])
    assert out["HSM_FLAG_PSF_2"][0] == 1  # pre-fill, never overwritten

    # But the identity of slot 2 is still recorded.
    assert out["EXP_ID_2"][0] == 2358123
    assert out["CCD_2"][0] == 21


def test_save_psf_data_fills_sentinel_for_absent_epochs(tmp_path):
    """Unused epoch slots and "empty" objects keep the -1 sentinel.

    ``max_epoch`` (from the largest geometric sexcat ``N_EPOCH``) can
    exceed a given object's own epoch count, and an object the PSF
    catalogue reports no epochs at all for is marked ``"empty"``; both
    cases must leave EXP_ID_n/CCD_n at -1, an exposure ID / CCD number no
    real epoch can have.
    """
    galaxy_psf_path = tmp_path / "galaxy_psf.sqlite"
    per_obj = {
        101: {"2113864-7": _psf_epoch(0.01, 0.02, 0.5)},
        303: "empty",
    }
    _write_galaxy_psf_cat(galaxy_psf_path, per_obj)

    # Geometric max(N_EPOCH) = 2 -> 3 slots; obj 101 only fills slot 1.
    out = _run_save_psf(galaxy_psf_path, [101, 303], n_overlap=[1, 2])

    assert out["EXP_ID_1"][0] == 2113864
    assert out["CCD_1"][0] == 7
    for col in ("EXP_ID_2", "CCD_2", "EXP_ID_3", "CCD_3"):
        assert out[col][0] == -1, col
    for col in ("EXP_ID_1", "CCD_1", "EXP_ID_2", "CCD_2", "EXP_ID_3", "CCD_3"):
        assert out[col][1] == -1, col


def test_save_psf_data_n_epoch_counts_only_psf_validated_epochs(tmp_path):
    """N_EPOCH counts the epochs with an interpolated PSF, not the overlaps.

    Object 101 overlaps three CCDs but one failed PSF-model validation, so
    the PSF catalogue holds two epochs for it and ngmix sees two; object 202
    has an epoch whose PSF shape fit failed, which still counts (its PSF
    exists); object 303 overlaps a CCD with no validated PSF at all and is
    ``"empty"``.
    """
    galaxy_psf_path = tmp_path / "galaxy_psf.sqlite"
    per_obj = {
        101: {
            "2603237-12": _psf_epoch(0.01, 0.02, 0.5),
            "2603241-12": _psf_epoch(0.03, 0.04, 0.6),
        },
        202: {
            "2603237-13": _psf_epoch(0.05, 0.06, 0.7),
            "2603241-13": _psf_epoch(-10.0, -10.0, 0.0, flag=5),
        },
        303: "empty",
    }
    _write_galaxy_psf_cat(galaxy_psf_path, per_obj)

    out = _run_save_psf(galaxy_psf_path, [101, 202, 303], n_overlap=[3, 2, 1])

    npt.assert_array_equal(out["N_EPOCH"], [2, 2, 0])
    assert out["EXP_ID_3"][0] == -1


def test_save_psf_data_without_epoch_slots_writes_only_n_epoch(tmp_path):
    """With per-epoch slots off, ``N_EPOCH`` is the only column written."""
    galaxy_psf_path = tmp_path / "galaxy_psf.sqlite"
    per_obj = {
        101: {"2603237-12": _psf_epoch(0.01, 0.02, 0.5)},
        303: "empty",
    }
    _write_galaxy_psf_cat(galaxy_psf_path, per_obj)

    out = _run_save_psf(
        galaxy_psf_path, [101, 303], n_overlap=[2, 1], epoch_slots=False
    )

    assert list(out) == ["N_EPOCH"]
    npt.assert_array_equal(out["N_EPOCH"], [1, 0])


# --- make_cat_runner: end-to-end catalogue assembly ---


def _write_sex_like_cat(path, data):
    """Write a minimal SExtractor-format FITS (data lives at HDU index 2)."""
    fits.HDUList(
        [
            fits.PrimaryHDU(),
            fits.BinTableHDU(name="LDAC_IMHEAD"),
            fits.BinTableHDU(data, name="LDAC_OBJECTS"),
        ]
    ).writeto(str(path), overwrite=True)


def _numbered_data(obj_ids):
    """A minimal SExtractor-like table: ``NUMBER`` plus one position field."""
    return np.array(
        [(oid, 0.0) for oid in obj_ids],
        dtype=[("NUMBER", "i8"), ("X_IMAGE", "f8")],
    )


def test_make_cat_runner_ships_every_detection_unclassified(tmp_path):
    """The runner assembles every detection, with no star/galaxy column.

    Star/galaxy separation happens downstream, so the final catalogue keeps
    each SExtractor object and carries no classification or spread-model
    column.
    """
    obj_ids = [1, 2, 3]
    tile_sexcat_path = tmp_path / "tile_sexcat-350-100.fits"
    _write_sex_like_cat(tile_sexcat_path, _numbered_data(obj_ids))
    galaxy_psf_path = tmp_path / "galaxy_psf-350-100.sqlite"
    _write_galaxy_psf_cat(galaxy_psf_path, dict.fromkeys(obj_ids, "empty"))
    ngmix_path = tmp_path / "ngmix-350-100.fits"
    _write_ngmix_cat(ngmix_path, obj_ids)

    config = CustomParser()
    config.read_dict({"MAKE_CAT_RUNNER": {"SHAPE_MEASUREMENT_TYPE": "ngmix"}})

    result = make_cat_runner(
        [str(tile_sexcat_path), str(galaxy_psf_path), str(ngmix_path)],
        {"output": str(tmp_path)},
        "-350-100",
        config,
        "MAKE_CAT_RUNNER",
        _NullLogger(),
    )

    assert result == (None, None)
    with h5py.File(tmp_path / "final_cat-350-100.hdf5", "r") as cat:
        npt.assert_array_equal(cat["NUMBER"][()], obj_ids)
        assert not [name for name in cat if "SPREAD" in name]
        # No MASK_EXT_PATHS, no mask columns.
        assert not [name for name in cat if name.startswith("MASK_")]


def test_make_cat_runner_n_epoch_excludes_unvalidated_psf_epochs(tmp_path):
    """The final N_EPOCH counts only PSF-validated epochs.

    The SExtractor ``N_EPOCH`` counts every exposure CCD covering an object.
    Object 1 lies on two CCDs of which one failed PSF validation, object 3
    on one such CCD only. Without SAVE_PSF_DATA the runner still writes
    ``N_EPOCH`` from the PSF catalogue, and no per-epoch slot column.
    """
    obj_ids = [1, 2, 3]
    sexcat = np.array(
        list(zip(obj_ids, [2, 3, 1])),
        dtype=[("NUMBER", "i8"), ("N_EPOCH", "i8")],
    )
    tile_sexcat_path = tmp_path / "tile_sexcat-350-100.fits"
    _write_sex_like_cat(tile_sexcat_path, sexcat)
    galaxy_psf_path = tmp_path / "galaxy_psf-350-100.sqlite"
    _write_galaxy_psf_cat(
        galaxy_psf_path,
        {
            1: {"2603241-12": _psf_epoch(0.01, 0.02, 0.5)},
            2: {
                "2603237-12": _psf_epoch(0.05, 0.06, 0.7),
                "2603241-12": _psf_epoch(0.07, 0.08, 0.9),
                "2603246-12": _psf_epoch(0.03, 0.04, 0.6),
            },
            3: "empty",
        },
    )
    ngmix_path = tmp_path / "ngmix-350-100.fits"
    _write_ngmix_cat(ngmix_path, obj_ids)

    config = CustomParser()
    config.read_dict({"MAKE_CAT_RUNNER": {"SHAPE_MEASUREMENT_TYPE": "ngmix"}})
    assert make_cat_runner(
        [str(tile_sexcat_path), str(galaxy_psf_path), str(ngmix_path)],
        {"output": str(tmp_path)},
        "-350-100",
        config,
        "MAKE_CAT_RUNNER",
        _NullLogger(),
    ) == (None, None)

    with h5py.File(tmp_path / "final_cat-350-100.hdf5", "r") as cat:
        npt.assert_array_equal(cat["N_EPOCH"][()], [1, 3, 0])
        assert "N_EPOCH_OVERLAP" not in cat
        assert "EXP_ID_1" not in cat


def test_make_cat_runner_writes_one_hdf5_dataset_per_column(tmp_path):
    """All save stages land in one hdf5 file, one lzf dataset per column.

    Runs every stage (ngmix, per-epoch PSF slots, one external mask band).
    Columns keep the order the stages add them, a vector column is a 2-D
    dataset with one row per object, and each column keeps the dtype its
    stage built, in native byte order (the SExtractor input is big-endian
    FITS).
    """
    healsparse = pytest.importorskip("healsparse")

    obj_ids = [1, 2, 3]
    ra = np.array([10.0, 10.1, 200.0])
    dec = np.array([20.0, 20.1, -40.0])
    flux_aper = np.arange(9, dtype=np.float32).reshape(3, 3)
    sexcat = np.array(
        list(zip(obj_ids, [1, 2, 0], ra, dec, flux_aper)),
        dtype=[
            ("NUMBER", "i8"),
            ("N_EPOCH", "i8"),
            ("XWIN_WORLD", "f8"),
            ("YWIN_WORLD", "f8"),
            ("FLUX_APER", "f4", (3,)),
        ],
    )
    tile_sexcat_path = tmp_path / "tile_sexcat-350-100.fits"
    _write_sex_like_cat(tile_sexcat_path, sexcat)
    galaxy_psf_path = tmp_path / "galaxy_psf-350-100.sqlite"
    _write_galaxy_psf_cat(
        galaxy_psf_path,
        {
            1: {"2113864-7": _psf_epoch(0.01, 0.02, 0.5)},
            2: {
                "2113864-9": _psf_epoch(0.05, 0.06, 0.7),
                "2358123-21": _psf_epoch(0.07, 0.08, 0.9),
            },
            3: "empty",
        },
    )
    ngmix_path = tmp_path / "ngmix-350-100.fits"
    _write_ngmix_cat(ngmix_path, obj_ids)
    mask = healsparse.HealSparseMap.make_empty(32, 4096, np.int16, sentinel=-1)
    mask.update_values_pos(
        ra[:2], dec[:2], np.full(2, 64, dtype=np.int16), lonlat=True
    )
    mask_path = tmp_path / "mask_r.hsp"
    mask.write(str(mask_path))

    config = CustomParser()
    config.read_dict({"MAKE_CAT_RUNNER": {
        "SHAPE_MEASUREMENT_TYPE": "ngmix",
        "SAVE_PSF_DATA": "True",
        "N_EPOCH_SLOTS": "3",
        "MASK_EXT_PATHS": f"r:{mask_path}",
    }})
    out_dir = tmp_path / "output"
    out_dir.mkdir()
    assert make_cat_runner(
        [str(tile_sexcat_path), str(galaxy_psf_path), str(ngmix_path)],
        {"output": str(out_dir)}, "-350-100", config,
        "MAKE_CAT_RUNNER", _NullLogger(),
    ) == (None, None)
    assert [p.name for p in out_dir.iterdir()] == ["final_cat-350-100.hdf5"]

    with h5py.File(out_dir / "final_cat-350-100.hdf5", "r") as cat:
        names = list(cat)
        assert names[:7] == [
            "NUMBER", "N_EPOCH", "XWIN_WORLD", "YWIN_WORLD",
            "FLUX_APER", "TILE_ID", "TILE_UNIQUE_ID",
        ]
        assert (
            names.index("N_EPOCH")
            < names.index("NGMIX_G1_NOSHEAR")
            < names.index("HSM_G1_PSF_1")
        )
        assert names[-1] == "MASK_r"
        for name in names:
            assert cat[name].shape[0] == len(obj_ids), name
            assert cat[name].compression == "lzf", name
            assert cat[name].dtype.byteorder in "=|<", name
        assert cat["FLUX_APER"].shape == (3, 3)
        npt.assert_array_equal(cat["FLUX_APER"][()], flux_aper)
        assert cat["NUMBER"].dtype == np.int64
        assert cat["TILE_UNIQUE_ID"].dtype == np.int64
        assert cat["HSM_FLAG_PSF_1"].dtype == np.int16
        assert cat["EXP_ID_2"].dtype == np.int32
        npt.assert_array_equal(cat["NUMBER"][()], obj_ids)
        npt.assert_array_equal(cat["N_EPOCH"][()], [1, 2, 0])
        npt.assert_allclose(cat["HSM_G1_PSF_1"][()], [0.01, 0.05, -10.0])
        npt.assert_array_equal(cat["EXP_ID_2"][()], [-1, 2358123, -1])
        npt.assert_array_equal(cat["MASK_r"][()], [64, 64, -1])


@pytest.mark.parametrize("shear", SHEAR_EXTS)
@pytest.mark.parametrize("component", [0, 1], ids=["g1", "g2"])
@pytest.mark.parametrize("nonfinite", [np.nan, np.inf, -np.inf])
def test_galaxy_cut_rejects_each_nonfinite_shear_component(
    tmp_path, shear, component, nonfinite,
):
    """A non-finite component in any metacal type cannot pass the galaxy cut."""
    result = _metacal_result(1)
    result[shear.lower()]["g"] = [0.1, 0.2]
    result[shear.lower()]["g"][component] = nonfinite
    out = _serialise_then_merge(tmp_path, [result], np.array([1]))

    assert out["NGMIX_MCAL_FLAGS"][0] == LM_FUNC_NOTFINITE
    assert out["NGMIX_MCAL_TYPES_FAIL"][0] == 1
    assert out[f"NGMIX_FLAGS_{shear}"][0] == LM_FUNC_NOTFINITE


def test_galaxy_cut_admits_only_measured_objects(tmp_path):
    """Contracts mcal-flags-zero-means-measured and failure-sentinel-cut-semantics.

    Consumer-side invariant through the real write path: synthetic metacal
    results -> ``compile_results`` / ``save_results`` -> ``_save_ngmix_data``
    -> sp_validation's galaxy cut ``MCAL_FLAGS == 0 & MCAL_TYPES_FAIL == 0``.
    One object per way a fit can fail to be a measurement, plus objects
    ngmix never fit. Only the two clean objects may pass, and every passing
    row must carry a fitted shape in all five metacal types.
    """
    nan, inf = float("nan"), float("inf")
    results = {oid: _metacal_result(oid) for oid in range(1, 9)}
    # 1, 8: clean.
    results[2]["1p"] = {"flags": 0x8, "nfev": 5}  # fitter reported failure
    del results[3]["noshear"]  # type absent from the result
    del results[4]["2m"]["flags"]  # no flags key, shape present
    results[5]["noshear"]["g"] = [nan, nan]  # flags 0, non-finite shear
    results[6]["1m"]["g"] = [inf, 0.1]  # flags 0, infinite shear
    del results[7]["2p"]["g"]  # flags 0, no shear at all
    # 101, 102: in the final catalogue but never fit by ngmix.
    cat_ids = np.array([101, 1, 2, 3, 102, 4, 5, 6, 7, 8])

    out = _serialise_then_merge(tmp_path, list(results.values()), cat_ids)

    mcal_flags = np.asarray(out["NGMIX_MCAL_FLAGS"]).astype(np.int64)
    types_fail = np.asarray(out["NGMIX_MCAL_TYPES_FAIL"]).astype(np.int64)
    passed = (mcal_flags == 0) & (types_fail == 0)

    assert set(cat_ids[passed]) == {1, 8}, (
        "mcal-flags-zero-means-measured / failure-sentinel-cut-semantics: the cut"
        f" admitted {sorted(set(cat_ids[passed]) - {1, 8})}"
    )
    assert np.all(np.asarray(out["NGMIX_N_EPOCH"])[passed] > 0)
    for shear in SHEAR_EXTS:
        for comp in ("G1", "G2"):
            shape = np.asarray(out[f"NGMIX_{comp}_{shear}"])[passed]
            assert np.all(np.isfinite(shape) & (shape != -10.0)), (
                f"NGMIX_{comp}_{shear}: a row passing the cut has no fitted shape"
            )

    # The three flag columns agree on every row, fitted or never fit:
    # MCAL_FLAGS is the OR and MCAL_TYPES_FAIL the count of FLAGS_<SHEAR>.
    type_flags = np.array(
        [np.asarray(out[f"NGMIX_FLAGS_{shear}"]) for shear in SHEAR_EXTS]
    ).astype(np.int64)
    npt.assert_array_equal(mcal_flags, np.bitwise_or.reduce(type_flags))
    npt.assert_array_equal(types_fail, np.count_nonzero(type_flags, axis=0))

    # Failures carry ngmix's own bits: the fitter's flags pass through, and
    # every no-finite-shear case, never-fit objects included, reads
    # LM_FUNC_NOTFINITE. No bit outside ngmix.flags is ever set.
    expected = {oid: LM_FUNC_NOTFINITE for oid in (3, 4, 5, 6, 7, 101, 102)}
    expected.update({1: 0, 8: 0, 2: 0x8})
    npt.assert_array_equal(mcal_flags, [expected[oid] for oid in cat_ids])
    ngmix_bits = functools.reduce(operator.or_, NAME_MAP)
    assert not np.any(type_flags & ~ngmix_bits)


# --- _save_psf_data: fixed per-epoch slot count (N_EPOCH_SLOTS) ---

# Per-family empty-slot sentinel: what a slot holds when no epoch fills it.
_PSF_SLOT_SENTINELS = {
    "HSM_G1_PSF": -10.0,
    "HSM_G2_PSF": -10.0,
    "HSM_T_PSF": 0.0,
    "HSM_FLAG_PSF": 1,
    "HSM_M4_1_PSF": -10.0,
    "HSM_M4_2_PSF": -10.0,
    "HSM_RHO4_PSF": -1.0,
    "EXP_ID": -1,
    "CCD": -1,
}


def _slot_numbers(out, family):
    """The slot numbers ``n`` present in ``out`` for ``<family>_n`` columns."""
    prefix = f"{family}_"
    return sorted(
        int(col[len(prefix):])
        for col in out
        if col.startswith(prefix) and col[len(prefix):].isdigit()
    )


def test_save_psf_data_fixed_slots_pad_every_family(tmp_path):
    """N_EPOCH_SLOTS fixes the slot count for every family, sentinel-padded.

    The campaign merge needs one schema across tiles, so a tile whose
    objects all have far fewer epochs than N_EPOCH_SLOTS still writes
    exactly slots 1..N_EPOCH_SLOTS for each per-epoch family, and every
    slot no epoch fills holds that family's own sentinel.
    """
    galaxy_psf_path = tmp_path / "galaxy_psf.sqlite"
    per_obj = {
        101: {"2113864-7": _psf_epoch(0.01, 0.02, 0.5)},
        202: {
            "2113864-9": _psf_epoch(0.05, 0.06, 0.7),
            "2358123-21": _psf_epoch(0.07, 0.08, 0.9),
        },
        303: "empty",
    }
    _write_galaxy_psf_cat(galaxy_psf_path, per_obj)

    # Driven through ``process``, the entry point the runner calls.
    n_slots = 7
    out = {
        "NUMBER": np.array([101, 202, 303]),
        "N_EPOCH": np.array([1, 2, 0]),
    }
    sc = SaveCatalogue(out, 3, _NullLogger())
    assert sc.process("psf", str(galaxy_psf_path), n_epoch_slots=n_slots) is None

    n_filled = [1, 2, 0]
    for family, sentinel in _PSF_SLOT_SENTINELS.items():
        assert _slot_numbers(out, family) == list(range(1, n_slots + 1)), family
        for row, filled in enumerate(n_filled):
            for n in range(filled + 1, n_slots + 1):
                assert out[f"{family}_{n}"][row] == sentinel, (family, row, n)


def test_save_psf_data_more_epochs_than_slots_raises(tmp_path):
    """An object with more epochs than N_EPOCH_SLOTS raises, never truncates.

    Dropping the extra epochs would silently lose per-epoch PSF data; the
    error names the object, its N_EPOCH and the slot count so the config
    can be fixed.
    """
    galaxy_psf_path = tmp_path / "galaxy_psf.sqlite"
    per_obj = {
        101: {"2113864-7": _psf_epoch(0.01, 0.02, 0.5)},
        404: {
            "2113864-7": _psf_epoch(0.01, 0.02, 0.5),
            "2229900-13": _psf_epoch(0.03, 0.04, 0.6),
            "2358123-21": _psf_epoch(0.07, 0.08, 0.9),
        },
    }
    _write_galaxy_psf_cat(galaxy_psf_path, per_obj)

    with pytest.raises(ValueError) as excinfo:
        _run_save_psf(
            galaxy_psf_path, [101, 404], n_overlap=[1, 3], n_epoch_slots=2
        )
    msg = str(excinfo.value)
    assert "404" in msg
    assert "N_EPOCH=3" in msg
    assert "N_EPOCH_SLOTS=2" in msg


def test_save_psf_data_exactly_slots_epochs_fits(tmp_path):
    """An object whose epochs exactly fill N_EPOCH_SLOTS is written, not raised."""
    galaxy_psf_path = tmp_path / "galaxy_psf.sqlite"
    per_obj = {
        404: {
            "2113864-7": _psf_epoch(0.01, 0.02, 0.5),
            "2229900-13": _psf_epoch(0.03, 0.04, 0.6),
        },
    }
    _write_galaxy_psf_cat(galaxy_psf_path, per_obj)

    out = _run_save_psf(galaxy_psf_path, [404], n_overlap=[2], n_epoch_slots=2)

    assert _slot_numbers(out, "EXP_ID") == [1, 2]
    assert out["EXP_ID_2"][0] == 2229900
    npt.assert_allclose(out["HSM_G1_PSF_2"], [0.03])


def test_save_psf_data_unset_slots_uses_tile_max_n_overlap_plus_one(tmp_path):
    """Without N_EPOCH_SLOTS use the sexcat's geometric max(N_EPOCH) + 1."""
    galaxy_psf_path = tmp_path / "galaxy_psf.sqlite"
    per_obj = {
        101: {"2113864-7": _psf_epoch(0.01, 0.02, 0.5)},
        202: {
            "2113864-9": _psf_epoch(0.05, 0.06, 0.7),
            "2229900-13": _psf_epoch(0.03, 0.04, 0.6),
            "2358123-21": _psf_epoch(0.07, 0.08, 0.9),
        },
    }
    _write_galaxy_psf_cat(galaxy_psf_path, per_obj)

    out = _run_save_psf(galaxy_psf_path, [101, 202], n_overlap=[1, 3])

    for family in _PSF_SLOT_SENTINELS:
        assert _slot_numbers(out, family) == [1, 2, 3, 4], family


def test_save_psf_data_fixed_slots_keep_epoch_alignment(tmp_path):
    """Under padding, slot n of every family still names the same epoch.

    Epochs fill slots 1..k in the galaxy_psf key order and padding sits
    only in slots k+1..N_EPOCH_SLOTS, for every family alike; a flagged
    epoch keeps its identity in its own slot while its HSM columns stay
    at the sentinel.
    """
    galaxy_psf_path = tmp_path / "galaxy_psf.sqlite"
    epochs = [
        # (key, g1, g2, t, flag) in deliberately unsorted key order
        ("2358123-21", 0.07, 0.08, 0.9, 0),
        ("2113864-9", -10.0, -10.0, 0.0, 5),
        ("2229900-13", 0.03, 0.04, 0.6, 0),
    ]
    per_obj = {
        505: {key: _psf_epoch(g1, g2, t, flag) for key, g1, g2, t, flag in epochs},
    }
    _write_galaxy_psf_cat(galaxy_psf_path, per_obj)

    n_slots = 6
    out = _run_save_psf(
        galaxy_psf_path, [505], n_overlap=[3], n_epoch_slots=n_slots
    )

    for n, (key, g1, g2, t, flag) in enumerate(epochs, start=1):
        exp_id, ccd = (int(part) for part in key.split("-"))
        assert out[f"EXP_ID_{n}"][0] == exp_id, n
        assert out[f"CCD_{n}"][0] == ccd, n
        if flag == 0:
            npt.assert_allclose(out[f"HSM_G1_PSF_{n}"], [g1])
            npt.assert_allclose(out[f"HSM_G2_PSF_{n}"], [g2])
            npt.assert_allclose(out[f"HSM_T_PSF_{n}"], [t])
            assert out[f"HSM_FLAG_PSF_{n}"][0] == 0, n
        else:
            npt.assert_allclose(out[f"HSM_G1_PSF_{n}"], [-10.0])
            assert out[f"HSM_FLAG_PSF_{n}"][0] == 1, n
    for n in range(len(epochs) + 1, n_slots + 1):
        for family, sentinel in _PSF_SLOT_SENTINELS.items():
            assert out[f"{family}_{n}"][0] == sentinel, (family, n)


def test_save_psf_data_carries_fourth_moments_per_epoch(tmp_path):
    """HSM_M4_1/M4_2/RHO4_PSF_n ride the same slots as HSM_G1_PSF_n.

    The fourth-moment columns psfex_interp writes into SHAPES (shapepipe#697)
    land per epoch; a SHAPES dict without them (MCCD, or an older producer)
    leaves the slot at its out-of-range fill, as does an unused slot.
    """
    galaxy_psf_path = tmp_path / "galaxy_psf.sqlite"
    per_obj = {
        101: {
            "2113864-7": _psf_epoch(0.01, 0.02, 0.5, m4=(0.11, -0.22, 2.05)),
            "2113865-3": _psf_epoch(0.03, 0.04, 0.6),
        },
    }
    _write_galaxy_psf_cat(galaxy_psf_path, per_obj)

    out = _run_save_psf(galaxy_psf_path, [101], n_overlap=[2])

    npt.assert_allclose(out["HSM_M4_1_PSF_1"], [0.11])
    npt.assert_allclose(out["HSM_M4_2_PSF_1"], [-0.22])
    npt.assert_allclose(out["HSM_RHO4_PSF_1"], [2.05])
    for n in (2, 3):
        npt.assert_allclose(out[f"HSM_M4_1_PSF_{n}"], [-10.0])
        npt.assert_allclose(out[f"HSM_M4_2_PSF_{n}"], [-10.0])
        npt.assert_allclose(out[f"HSM_RHO4_PSF_{n}"], [-1.0])
    npt.assert_allclose(out["HSM_G1_PSF_2"], [0.03])


def _write_sexcat(path, number):
    """Write a minimal FITS-LDAC tile catalogue (``LDAC_OBJECTS`` at HDU 2)."""
    n_obj = len(number)
    cols = [
        fits.Column(name="NUMBER", format="J", array=np.asarray(number)),
        fits.Column(name="XWIN_WORLD", format="D", array=np.zeros(n_obj)),
        fits.Column(
            name="VIGNET", format="4E", dim="(2,2)",
            array=np.zeros((n_obj, 2, 2), dtype=np.float32),
        ),
    ]
    fits.HDUList([
        fits.PrimaryHDU(),
        fits.BinTableHDU.from_columns(
            [fits.Column(name="Field Header Card", format="80A",
                         array=np.array([""]))],
            name="LDAC_IMHEAD",
        ),
        fits.BinTableHDU.from_columns(cols, name="LDAC_OBJECTS"),
    ]).writeto(path, overwrite=True)


def test_read_sextractor_data_adds_tile_unique_id(tmp_path):
    """SExtractor-mode final catalogue carries tile_id * 10**6 + NUMBER.

    ``NUMBER`` is deliberately gapped and unsorted: the ID is built from the
    value each row carries, not from its position.
    """
    number = np.array([7, 3, 12, 999999])
    sexcat = tmp_path / "sexcat-301-279.fits"
    _write_sexcat(sexcat, number)

    data = make_cat.read_sextractor_data(str(sexcat))

    assert list(data) == ["NUMBER", "XWIN_WORLD", "TILE_ID", "TILE_UNIQUE_ID"]
    npt.assert_array_equal(data["NUMBER"], number)
    assert data["TILE_UNIQUE_ID"].dtype.newbyteorder("=") == np.int64
    npt.assert_array_equal(data["TILE_UNIQUE_ID"], 301279 * 10**6 + number)
    npt.assert_allclose(data["TILE_ID"], 301.279)


def test_read_sextractor_data_refuses_number_beyond_id_range(tmp_path):
    """A NUMBER that would overflow into the tile digits raises."""
    sexcat = tmp_path / "sexcat-301-279.fits"
    _write_sexcat(sexcat, np.array([1, 10**6]))

    with pytest.raises(cfis.CfisError):
        make_cat.read_sextractor_data(str(sexcat))
