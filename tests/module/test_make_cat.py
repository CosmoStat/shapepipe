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
    the sentinels: 0 for sizes/flags, -10 for ellipticities, 1e30 for
    ``T_ERR``, -1 for flux/mag errors.
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


def test_save_ngmix_data_matches_module_serialised_catalogue(tmp_path):
    """End-to-end: make_cat reads a catalogue ngmix itself serialised.

    Rather than hand-rolling the FITS, drive ngmix's own
    ``compile_results`` + ``save_results`` (the real write path), then read it
    back through ``_save_ngmix_data``. This guards the make_cat reader against
    drift in the ngmix key set: the keys must line up end to end.
    """
    obj_ids = [5, 7]
    results = []
    for oid in obj_ids:
        row = _ngmix_row(oid)
        per_type = {
            "nfev": row["nfev_fit"],
            "g": [row["g1"], row["g2"]],
            "g_cov": np.diag([row["g1_err"] ** 2, row["g2_err"] ** 2]),
            "T": row["T"], "T_err": row["T_err"],
            "flux": row["flux"], "flux_err": row["flux_err"],
            "s2n": row["s2n"], "flags": 0,
        }
        res = {
            "obj_id": oid,
            "n_epoch_model": row["n_epoch_model"],
            "mcal_types_fail": row["mcal_types_fail"],
            "neighbour_flag": row["neighbour_flag"],
            "mcal_flags": row["mcal_flags"],
        }
        for key in NGMIX_KEYS:
            if key.endswith("_psf_orig") or key.endswith("_psf_reconv"):
                res[key] = row[key]
        res.update({name: dict(per_type) for name in SHEAR_EXTS_LOWER})
        results.append(res)

    ngmix_inst = object.__new__(Ngmix)
    ngmix_inst._zero_point = 30.0
    ngmix_inst._output_dir = str(tmp_path)
    ngmix_inst._file_number_string = "-0"
    out_dict = ngmix_inst.compile_results(results)
    ngmix_inst.save_results(out_dict)

    ngmix_path = ngmix_inst.get_output_path(str(tmp_path))
    out = _run_save_ngmix(ngmix_path, obj_ids)

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


def _run_save_psf(galaxy_psf_path, obj_id, n_epoch, n_epoch_slots=None):
    """Drive ``_save_psf_data`` and return its populated output dict."""
    inst = object.__new__(SaveCatalogue)
    inst._obj_id = np.asarray(obj_id)
    inst._output_dict = {}
    inst._final_cat_file = _FinalCatStub(n_epoch)

    inst._save_psf_data(str(galaxy_psf_path), n_epoch_slots=n_epoch_slots)
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


# --- _save_psf_data: fixed per-epoch slot count (N_EPOCH_SLOTS) ---

# Per-family empty-slot sentinel: what a slot holds when no epoch fills it.
_PSF_SLOT_SENTINELS = {
    "HSM_G1_PSF": -10.0,
    "HSM_G2_PSF": -10.0,
    "HSM_T_PSF": 0.0,
    "HSM_FLAG_PSF": 1,
    "EXP_ID": -1,
    "CCD": -1,
}


class _ProcessCatStub(_FinalCatStub):
    """FITSCatalogue stand-in for ``SaveCatalogue.process``; records add_col."""

    def __init__(self, obj_id, n_epoch):
        super().__init__(n_epoch)
        self._number = np.asarray(obj_id)
        self.cols = {}

    def open(self):
        pass

    def close(self):
        pass

    def get_data(self):
        return {"NUMBER": self._number, "N_EPOCH": self._n_epoch}

    def add_col(self, name, data):
        self.cols[name] = data


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
    cat = _ProcessCatStub([101, 202, 303], n_epoch=[1, 2, 0])
    sc = SaveCatalogue(cat, 3, _NullLogger())
    assert sc.process("psf", str(galaxy_psf_path), n_epoch_slots=n_slots) is None
    out = cat.cols

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
            galaxy_psf_path, [101, 404], n_epoch=[1, 3], n_epoch_slots=2
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

    out = _run_save_psf(galaxy_psf_path, [404], n_epoch=[2], n_epoch_slots=2)

    assert _slot_numbers(out, "EXP_ID") == [1, 2]
    assert out["EXP_ID_2"][0] == 2229900
    npt.assert_allclose(out["HSM_G1_PSF_2"], [0.03])


def test_save_psf_data_unset_slots_uses_tile_max_n_epoch_plus_one(tmp_path):
    """Without N_EPOCH_SLOTS the slot count is the tile's max(N_EPOCH) + 1."""
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

    out = _run_save_psf(galaxy_psf_path, [101, 202], n_epoch=[1, 3])

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
        galaxy_psf_path, [505], n_epoch=[3], n_epoch_slots=n_slots
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
