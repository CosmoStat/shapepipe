"""UNIT TESTS FOR MODULE PACKAGE: MAKE_CAT.

Drives ``SaveCatalogue._save_ngmix_data`` against a synthetic ngmix catalogue
to lock in the column grammar ``ESTIMATOR_COMPONENT[_ERR]_OBJECT[_SHEAR]``
(shapepipe#749, #761): galaxy is the implicit default object and carries no
``GAL`` token (``NGMIX_G1_NOSHEAR``, never ``NGMIX_G1_GAL_NOSHEAR``), while
the PSF families keep an explicit object token —
``NGMIX_<COMPONENT>[_ERR]_<OBJECT>_<SHEAR>`` for ``PSF_ORIG``/``PSF_RECONV`` —
plus the three OBJECT/SHEAR-less metadata columns (``NGMIX[m]_MCAL_FLAGS``,
``NGMIX_N_EPOCH``, ``NGMIX_MCAL_TYPES_FAIL``). The original image PSF
(``PSF_ORIG``) and the metacal reconvolution kernel (``PSF_RECONV``) are
independent fits of *different* PSFs, no longer the single aliased value of
the pre-#749 code.
"""

import numpy as np
import numpy.testing as npt
from astropy.io import fits

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
    for col in ("NGMIX_MCAL_FLAGS", "NGMIX_N_EPOCH", "NGMIX_MCAL_TYPES_FAIL"):
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
