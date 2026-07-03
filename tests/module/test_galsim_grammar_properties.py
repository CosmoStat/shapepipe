"""PROPERTY-BASED TESTS FOR THE GALSIM GALAXY/PSF COLUMN GRAMMAR.

Drives ``SaveCatalogue._save_galsim_shapes`` (the make_cat writer) to pin the
renamed galsim grammar (shapepipe#761): the packed length-2
``GALSIM_GAL_ELL_*`` / ``GALSIM_PSF_ELL_*`` columns are split into scalar
``GALSIM_G1_*`` / ``GALSIM_G2_*`` (g-type — galsim is correctly g-type, no
e/g value fix needed here), the galaxy carries no ``GAL`` token (implicit
default object, same rule as ngmix), and the raw ``GALSIM_GAL_SIGMA_*`` /
``GALSIM_PSF_SIGMA_*`` columns are replaced by the single stored size
``GALSIM_T_*`` (``T = 2 sigma^2``, via ``cs_util.size.sigma_to_T``).

The ``ORIGINAL_PSF`` extension is a deliberate asymmetry, preserved rather
than normalized away: its PSF columns are sourced from the *galaxy's*
uncorrected fields (``gal_uncorr_g1/g2``, ``gal_sigma``), not the
per-extension ``psf_g1/g2``/``psf_sigma`` fields every other extension uses.
"""

import re
import tempfile
from pathlib import Path

import numpy as np
import numpy.testing as npt
from astropy.io import fits
from cs_util import size as cs_size
from hypothesis import given, settings
from hypothesis import strategies as st

from shapepipe.modules.make_cat_package.make_cat import SaveCatalogue

# ---------------------------------------------------------------------------
# Harness
# ---------------------------------------------------------------------------

# Source columns _save_galsim_shapes reads per extension. The ORIGINAL_PSF
# extension only consumes the gal_uncorr_g1/g2 + gal_sigma subset, but the
# real galsim catalogue carries the full set on every extension, so the
# fixture mirrors that rather than special-casing the row shape.
GALSIM_KEYS = [
    "id",
    "gal_g1", "gal_g2", "gal_g1_err", "gal_g2_err",
    "gal_uncorr_g1", "gal_uncorr_g2", "gal_sigma",
    "psf_g1", "psf_g2", "psf_sigma",
    "gal_flux", "gal_flux_err", "gal_mag", "gal_mag_err",
    "gal_flag", "gal_resolution",
]
INT_KEYS = {"id", "gal_flag"}

# Per-key multiplier so a single scalar ``seed`` fans out into twelve
# distinct values — a swapped source (e.g. psf_g1 <-> gal_uncorr_g1) breaks
# at least one equality below.
_MULT = {
    "gal_g1": 1.0, "gal_g2": 1.3, "gal_g1_err": 1.7, "gal_g2_err": 2.1,
    "gal_uncorr_g1": 2.6, "gal_uncorr_g2": 3.1, "gal_sigma": 0.4,
    "psf_g1": 0.6, "psf_g2": 0.8, "psf_sigma": 0.5,
    "gal_flux": 10.0, "gal_flux_err": 1.0, "gal_mag": 5.0,
    "gal_mag_err": 0.2,
}


def _galsim_row(obj_id, seed):
    """One object's per-key values for one extension, scaled by ``seed``."""
    row = {"id": obj_id, "gal_flag": 0, "gal_resolution": 0.5}
    row.update({key: mult * seed for key, mult in _MULT.items()})
    return row


def _write_galsim_cat(path, obj_ids, seed_by_ext):
    """Write a synthetic galsim FITS, one extension per ``seed_by_ext`` key.

    Extension order is preserved (dict order), matching how
    ``get_ext_name()[1:]`` drives ``self._key_ends`` in production.
    """
    hdus = [fits.PrimaryHDU()]
    for ext, seed in seed_by_ext.items():
        rows = [_galsim_row(oid, seed) for oid in obj_ids]
        cols = [
            fits.Column(
                name=key,
                format="K" if key in INT_KEYS else "D",
                array=np.array([row[key] for row in rows]),
            )
            for key in GALSIM_KEYS
        ]
        hdus.append(fits.BinTableHDU.from_columns(cols, name=ext))
    fits.HDUList(hdus).writeto(path, overwrite=True)


def _run_save_galsim(galsim_path, obj_ids):
    """Drive ``_save_galsim_shapes`` and return its populated output dict."""
    inst = object.__new__(SaveCatalogue)
    inst._obj_id = np.asarray(obj_ids)
    inst._output_dict = {}
    inst._save_galsim_shapes(str(galsim_path))
    return inst._output_dict


# ---------------------------------------------------------------------------
# Strategies
# ---------------------------------------------------------------------------

obj_ids_strategy = st.lists(
    st.integers(min_value=1, max_value=10_000), min_size=1, max_size=5,
    unique=True,
)
# A finite seed bounded away from zero, so every per-key value it scales is
# itself non-zero and the family-vs-family / key-vs-key comparisons below are
# never vacuously satisfied by an accidental 0 == 0.
_seed = st.floats(
    min_value=0.05, max_value=5.0, allow_nan=False, allow_infinity=False,
)


# ---------------------------------------------------------------------------
# (a) Scalar split + size routing, regular extension
# ---------------------------------------------------------------------------

@settings(deadline=None, max_examples=25)
@given(obj_ids=obj_ids_strategy, seed=_seed)
def test_regular_extension_scalar_split_and_size(obj_ids, seed):
    """A non-``ORIGINAL_PSF`` extension routes each source to its own column.

    ``GALSIM_G1/G2_<K>`` from ``gal_g1/g2``, the ``_ERR``/``_UNCORR``
    variants from their namesake sources, ``GALSIM_T_<K>`` from
    ``sigma_to_T(gal_sigma)``, and the PSF family
    (``GALSIM_G1/G2_PSF_<K>``, ``GALSIM_T_PSF_<K>``) from ``psf_g1/g2`` /
    ``sigma_to_T(psf_sigma)`` — never from the galaxy's own fields.
    """
    ext = "NOSHEAR"
    path = tempfile.NamedTemporaryFile(suffix=".fits", delete=False).name
    _write_galsim_cat(path, obj_ids, {ext: seed})
    out = _run_save_galsim(path, obj_ids)
    Path(path).unlink()

    n = len(obj_ids)
    row = _galsim_row(obj_ids[0], seed)
    npt.assert_allclose(out[f"GALSIM_G1_{ext}"], [row["gal_g1"]] * n)
    npt.assert_allclose(out[f"GALSIM_G2_{ext}"], [row["gal_g2"]] * n)
    npt.assert_allclose(out[f"GALSIM_G1_ERR_{ext}"], [row["gal_g1_err"]] * n)
    npt.assert_allclose(out[f"GALSIM_G2_ERR_{ext}"], [row["gal_g2_err"]] * n)
    npt.assert_allclose(
        out[f"GALSIM_G1_UNCORR_{ext}"], [row["gal_uncorr_g1"]] * n
    )
    npt.assert_allclose(
        out[f"GALSIM_G2_UNCORR_{ext}"], [row["gal_uncorr_g2"]] * n
    )
    npt.assert_allclose(
        out[f"GALSIM_T_{ext}"], [cs_size.sigma_to_T(row["gal_sigma"])] * n
    )
    npt.assert_allclose(out[f"GALSIM_G1_PSF_{ext}"], [row["psf_g1"]] * n)
    npt.assert_allclose(out[f"GALSIM_G2_PSF_{ext}"], [row["psf_g2"]] * n)
    npt.assert_allclose(
        out[f"GALSIM_T_PSF_{ext}"], [cs_size.sigma_to_T(row["psf_sigma"])] * n
    )
    # The galaxy and PSF families must stay distinct (un-aliased), and the
    # raw-sigma value must never leak through unconverted.
    assert not np.allclose(out[f"GALSIM_T_{ext}"], row["gal_sigma"])
    assert not np.allclose(out[f"GALSIM_G1_{ext}"], out[f"GALSIM_G1_PSF_{ext}"])

    npt.assert_allclose(out[f"GALSIM_FLUX_{ext}"], [row["gal_flux"]] * n)
    npt.assert_allclose(out[f"GALSIM_FLUX_ERR_{ext}"], [row["gal_flux_err"]] * n)
    npt.assert_allclose(out[f"GALSIM_MAG_{ext}"], [row["gal_mag"]] * n)
    npt.assert_allclose(out[f"GALSIM_MAG_ERR_{ext}"], [row["gal_mag_err"]] * n)
    npt.assert_allclose(out[f"GALSIM_RES_{ext}"], [row["gal_resolution"]] * n)


# ---------------------------------------------------------------------------
# (b) ORIGINAL_PSF asymmetry: sourced from the galaxy's uncorrected fields
# ---------------------------------------------------------------------------

@settings(deadline=None, max_examples=25)
@given(obj_ids=obj_ids_strategy, seed_main=_seed, seed_origpsf=_seed)
def test_original_psf_extension_sourced_from_gal_uncorr(
    obj_ids, seed_main, seed_origpsf
):
    """``ORIGINAL_PSF`` columns trace to ``gal_uncorr_*``/``gal_sigma``, NOT
    ``psf_*`` — the one preserved asymmetry in the galsim grammar.

    The two extensions are written with independent seeds, and the
    ``ORIGINAL_PSF`` row's own ``psf_g1/g2``/``psf_sigma`` are deliberately
    distinct from its ``gal_uncorr_*``/``gal_sigma`` (the seed multipliers
    differ), so a regression that wired the regular per-extension ``psf_*``
    sourcing into ``ORIGINAL_PSF`` would be caught, not masked by equal
    values.
    """
    seed_by_ext = {"NOSHEAR": seed_main, "ORIGINAL_PSF": seed_origpsf}
    path = tempfile.NamedTemporaryFile(suffix=".fits", delete=False).name
    _write_galsim_cat(path, obj_ids, seed_by_ext)
    out = _run_save_galsim(path, obj_ids)
    Path(path).unlink()

    n = len(obj_ids)
    op_row = _galsim_row(obj_ids[0], seed_origpsf)
    npt.assert_allclose(
        out["GALSIM_G1_PSF_ORIGINAL_PSF"], [op_row["gal_uncorr_g1"]] * n
    )
    npt.assert_allclose(
        out["GALSIM_G2_PSF_ORIGINAL_PSF"], [op_row["gal_uncorr_g2"]] * n
    )
    npt.assert_allclose(
        out["GALSIM_T_PSF_ORIGINAL_PSF"],
        [cs_size.sigma_to_T(op_row["gal_sigma"])] * n,
    )
    # Not sourced from this extension's own psf_g1/g2/psf_sigma (those are
    # scaled by the same seed but through different multipliers, so equality
    # here would mean the asymmetry was lost).
    assert not np.allclose(
        out["GALSIM_G1_PSF_ORIGINAL_PSF"], [op_row["psf_g1"]] * n
    )
    assert not np.allclose(
        out["GALSIM_T_PSF_ORIGINAL_PSF"],
        [cs_size.sigma_to_T(op_row["psf_sigma"])] * n,
    )


# ---------------------------------------------------------------------------
# (c) Every emitted column matches the new grammar; no legacy tokens survive
# ---------------------------------------------------------------------------

@settings(deadline=None, max_examples=25)
@given(obj_ids=obj_ids_strategy, seed=_seed)
def test_emitted_columns_match_grammar_no_legacy_tokens(obj_ids, seed):
    """Every ``GALSIM_*`` column matches the scalar grammar; nothing packed.

    No ``GALSIM_GAL_ELL*``/``GALSIM_PSF_ELL*`` (packed length-2 ellipticity)
    and no ``GALSIM_*_SIGMA*`` (raw size) survives — both replaced by scalar
    ``G1``/``G2`` components and the single stored size ``T``. The galaxy
    also carries no ``GAL`` token (implicit default object, same rule as
    ngmix).
    """
    ext_names = ["NOSHEAR", "ORIGINAL_PSF"]
    seed_by_ext = {ext: seed * (i + 1) for i, ext in enumerate(ext_names)}
    path = tempfile.NamedTemporaryFile(suffix=".fits", delete=False).name
    _write_galsim_cat(path, obj_ids, seed_by_ext)
    out = _run_save_galsim(path, obj_ids)
    Path(path).unlink()

    assert out, "writer produced no columns"

    _component = "G1|G2|T"
    _qual = "(?:_ERR|_UNCORR)?"
    _ext = "|".join(re.escape(e) for e in ext_names)
    grammar_re = re.compile(
        rf"^GALSIM_(?:{_component}){_qual}(?:_PSF)?_(?:{_ext})$"
        rf"|^GALSIM_(?:FLUX|FLUX_ERR|MAG|MAG_ERR|FLAGS|RES)_(?:{_ext})$"
    )
    for col in out:
        assert grammar_re.match(col), f"off-grammar column emitted: {col!r}"
        assert "_ELL" not in col, col
        assert "SIGMA" not in col, col
        assert "_GAL_" not in col and not col.endswith("_GAL"), col
