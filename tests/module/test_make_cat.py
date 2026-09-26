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

import numpy as np
import numpy.testing as npt
import pytest
from astropy.io import fits
from sqlitedict import SqliteDict

from shapepipe.modules.make_cat_package.make_cat import SaveCatalogue
from shapepipe.modules.ngmix_package.ngmix import Ngmix


class _NullLogger:
    def info(self, *_args, **_kwargs):
        pass


# Per-object key set that compile_results emits into every shear-type
# extension (mirrors ngmix.compile_results ``names2``). Sentinel values are
# distinct per key so a mis-routed column is caught by value, and the orig vs
# reconv PSF families carry *different* sentinels so any aliasing surfaces.
NGMIX_KEYS = [
    "id",
    "n_epoch_model",
    "mcal_types_fail",
    "neighbour_flag",
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
        "mcal_types_fail": 0,
        "neighbour_flag": 1,
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
                format="K" if key in ("id", "n_epoch_model", "mcal_types_fail",
                                       "nfev_fit", "flags", "mcal_flags") else "D",
                array=np.array([row[key] for row in rows]),
            )
            for key in NGMIX_KEYS
        ]
        hdus.append(fits.BinTableHDU.from_columns(cols, name=ext))
    fits.HDUList(hdus).writeto(path, overwrite=True)


def _run_save_ngmix(ngmix_path, obj_id, cat_size_target=None):
    """Drive ``_save_ngmix_data`` and return its populated output dict."""
    inst = object.__new__(SaveCatalogue)
    inst._obj_id = np.asarray(obj_id)
    inst._output_dict = {}
    inst._cat_size_target = (
        len(inst._obj_id) if cat_size_target is None else cat_size_target
    )
    inst._w_log = _NullLogger()

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
        "NGMIX_NEIGHBOUR_FLAG",
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
        "neighbour_flag": row["neighbour_flag"],
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


def _psf_epoch(g1, g2, t, flag=0):
    """One epoch's interpolated-PSF HSM shape entry (psfex_interp's SHAPES dict)."""
    return {
        "SHAPES": {
            "HSM_G1_PSF": g1,
            "HSM_G2_PSF": g2,
            "HSM_T_PSF": t,
            "HSM_FLAG_PSF": flag,
        },
    }


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


class _FinalCatStub:
    """Stand-in for the FITSCatalogue ``_save_psf_data`` reads N_EPOCH from."""

    def __init__(self, n_epoch):
        self._n_epoch = np.asarray(n_epoch)

    def get_data(self):
        return {"N_EPOCH": self._n_epoch}


def _run_save_psf(galaxy_psf_path, obj_id, n_epoch):
    """Drive ``_save_psf_data`` and return its populated output dict."""
    inst = object.__new__(SaveCatalogue)
    inst._obj_id = np.asarray(obj_id)
    inst._output_dict = {}
    inst._final_cat_file = _FinalCatStub(n_epoch)

    inst._save_psf_data(str(galaxy_psf_path))
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

    out = _run_save_psf(galaxy_psf_path, [101], n_epoch=[2])

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

    out = _run_save_psf(galaxy_psf_path, [202], n_epoch=[2])

    npt.assert_allclose(out["HSM_G1_PSF_1"], [0.05])
    # Sentinel: the epoch-2 HSM fit failed, so the pre-fill value stands.
    npt.assert_allclose(out["HSM_G1_PSF_2"], [-10.0])
    assert out["HSM_FLAG_PSF_2"][0] == 1  # pre-fill, never overwritten

    # But the identity of slot 2 is still recorded.
    assert out["EXP_ID_2"][0] == 2358123
    assert out["CCD_2"][0] == 21


def test_save_psf_data_fills_sentinel_for_absent_epochs(tmp_path):
    """Unused epoch slots and "empty" objects keep the -1 sentinel.

    ``max_epoch`` (from the largest ``N_EPOCH`` across the catalogue) can
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

    # max(N_EPOCH) = 2 -> 3 slot columns, though obj 101 only fills slot 1.
    out = _run_save_psf(galaxy_psf_path, [101, 303], n_epoch=[1, 2])

    assert out["EXP_ID_1"][0] == 2113864
    assert out["CCD_1"][0] == 7
    for col in ("EXP_ID_2", "CCD_2", "EXP_ID_3", "CCD_3"):
        assert out[col][0] == -1, col
    for col in ("EXP_ID_1", "CCD_1", "EXP_ID_2", "CCD_2", "EXP_ID_3", "CCD_3"):
        assert out[col][1] == -1, col


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

    assert out["NGMIX_MCAL_FLAGS"][0] != 0
    assert out["NGMIX_MCAL_TYPES_FAIL"][0] == 1
    assert out[f"NGMIX_FLAGS_{shear}"][0] != 0


def test_galaxy_cut_admits_only_measured_objects(tmp_path):
    """Contracts mcal-flags-zero-means-measured and never-fit-is-not-clean.

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
        "mcal-flags-zero-means-measured / never-fit-is-not-clean: the cut"
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
