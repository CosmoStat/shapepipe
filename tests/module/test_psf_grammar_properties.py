"""PROPERTY-BASED TESTS FOR THE NGMIX PSF COLUMN GRAMMAR.

Drives ``SaveCatalogue._save_ngmix_data`` (the make_cat writer) under
hypothesis to pin three invariants of the renamed grammar
``NGMIX[m]_<COMPONENT>[_ERR][_<OBJECT>]_<SHEAR>`` (shapepipe#749, #761): the
galaxy is the implicit default object and carries NO ``OBJECT`` token
(``NGMIX_G1_NOSHEAR``, never ``NGMIX_G1_GAL_NOSHEAR``), plus the four
OBJECT/SHEAR-less metadata columns (``NGMIX[m]_MCAL_FLAGS``, ``NGMIX_N_EPOCH``,
``NGMIX_MCAL_TYPES_FAIL``, ``NGMIX_NEIGHBOUR_FLAG``), where the ORIGINAL image
PSF (``PSF_ORIG``) and the
metacal RECONVOLUTION kernel (``PSF_RECONV``) are independent fits of
*different* PSFs:

(a) un-aliasing holds for ALL inputs: for arbitrary distinct values fed as the
    orig vs reconv ngmix source columns, every emitted ``*_PSF_ORIG_*`` column
    traces to its orig source and every ``*_PSF_RECONV_*`` to its reconv
    source, and the two families are never crossed;
(b) every column name the writer emits matches the grammar regex (the OBJECT
    group is optional, so a regressed ``GAL`` token is rejected, not
    absorbed);
(c) every ``NGMIX_*`` token the shipped param file
    (``example/cfis/final_cat.param``) names is a column the writer can
    produce — writer/param-file consistency;
(d) the FULL frozen grammar (shapepipe#761) — ``ESTIMATOR_COMPONENT[_ERR]_
    OBJECT[_metacaltype]`` — holds across both estimator families: both
    families split ellipticity into ``G1``/``G2`` (g-type, HSM included since
    the return to reduced shear); every family stores exactly one size, ``T``; HSM
    keeps the singular ``FLAG`` token, ngmix keeps plural ``FLAGS``.
    This part is a static example-based check (not writer-driven), so it
    covers HSM without depending on that producer module — see
    ``test_psfex_interp.py`` for the writer-driven coverage of that family.

Setup for (a)-(c) mirrors ``tests/module/test_make_cat.py``: a synthetic
ngmix FITS with the five shear-type extensions carrying the exact
``compile_results`` key set, read back through ``_save_ngmix_data``.
"""

import re
from pathlib import Path

import numpy as np
import numpy.testing as npt
import pytest
from astropy.io import fits
from hypothesis import given, settings
from hypothesis import strategies as st

from shapepipe.modules.make_cat_package.make_cat import SaveCatalogue


# ---------------------------------------------------------------------------
# Harness (mirrors tests/module/test_make_cat.py)
# ---------------------------------------------------------------------------

class _NullLogger:
    def info(self, *_args, **_kwargs):
        pass


SHEAR_EXTS = ["1M", "1P", "2M", "2P", "NOSHEAR"]

# The full per-object key set ngmix.compile_results emits into every
# shear-type extension. Integer-typed keys carry FITS format "K", the rest
# "D"; everything else mirrors test_make_cat._ngmix_row.
INT_KEYS = {
    "id", "n_epoch_model", "mcal_types_fail", "neighbour_flag", "nfev_fit",
    "flags", "mcal_flags",
}
NGMIX_KEYS = [
    "id", "n_epoch_model", "mcal_types_fail", "neighbour_flag", "nfev_fit",
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

# The six orig-PSF source keys and the six reconv-PSF source keys, paired by
# the output column they feed. Each pair is (output-column infix, source key);
# the infix slots into f"NGMIX_{infix}_{shear}". Holding this mapping
# explicitly is what makes the un-aliasing property check each column against
# its OWN source rather than against a single family-wide value, so a within-
# family swap (e.g. G1 <-> T, or value <-> error) is caught too.
ORIG_COLS = {
    "G1_PSF_ORIG": "g1_psf_orig",
    "G2_PSF_ORIG": "g2_psf_orig",
    "G1_ERR_PSF_ORIG": "g1_err_psf_orig",
    "G2_ERR_PSF_ORIG": "g2_err_psf_orig",
    "T_PSF_ORIG": "T_psf_orig",
    "T_ERR_PSF_ORIG": "T_err_psf_orig",
}
RECONV_COLS = {
    "G1_PSF_RECONV": "g1_psf_reconv",
    "G2_PSF_RECONV": "g2_psf_reconv",
    "G1_ERR_PSF_RECONV": "g1_err_psf_reconv",
    "G2_ERR_PSF_RECONV": "g2_err_psf_reconv",
    "T_PSF_RECONV": "T_psf_reconv",
    "T_ERR_PSF_RECONV": "T_err_psf_reconv",
}


def _base_row(obj_id):
    """One object's per-key values; PSF keys overwritten by the caller."""
    return {
        "id": obj_id,
        "n_epoch_model": 3, "mcal_types_fail": 0, "neighbour_flag": 1,
        "nfev_fit": 7,
        "g1": 0.10, "g1_err": 0.011, "g2": -0.20, "g2_err": 0.022,
        "T": 0.30, "T_err": 0.033,
        "flux": 100.0, "flux_err": 1.0, "s2n": 55.0,
        "mag": 21.0, "mag_err": 0.05,
        "flags": 0, "mcal_flags": 0,
        "g1_psf_orig": 0.0041, "g2_psf_orig": -0.0031,
        "g1_err_psf_orig": 2e-5, "g2_err_psf_orig": 3e-5,
        "T_psf_orig": 0.071, "T_err_psf_orig": 8e-4,
        "g1_psf_reconv": 0.0012, "g2_psf_reconv": 0.0022,
        "g1_err_psf_reconv": 1e-5, "g2_err_psf_reconv": 1e-5,
        "T_psf_reconv": 0.092, "T_err_psf_reconv": 1e-3,
    }


def _write_ngmix_cat(path, rows):
    """Write a synthetic ngmix FITS with the five shear-type extensions.

    Each extension carries the exact compile_results key set; PSF families are
    object-level, so the same per-object values are written across extensions.
    """
    hdus = [fits.PrimaryHDU()]
    for ext in SHEAR_EXTS:
        cols = [
            fits.Column(
                name=key,
                format="K" if key in INT_KEYS else "D",
                array=np.array([row[key] for row in rows]),
            )
            for key in NGMIX_KEYS
        ]
        hdus.append(fits.BinTableHDU.from_columns(cols, name=ext))
    fits.HDUList(hdus).writeto(path, overwrite=True)


def _run_save_ngmix(ngmix_path, obj_ids, cat_size_target=None):
    """Drive _save_ngmix_data and return its populated output dict."""
    inst = object.__new__(SaveCatalogue)
    inst._obj_id = np.asarray(obj_ids)
    inst._output_dict = {}
    inst._cat_size_target = (
        len(inst._obj_id) if cat_size_target is None else cat_size_target
    )
    inst._w_log = _NullLogger()

    err_msg = inst._save_ngmix_data(str(ngmix_path))
    assert err_msg is None
    return inst._output_dict


# ---------------------------------------------------------------------------
# Grammar regex + producible-column model
# ---------------------------------------------------------------------------

# NGMIX[m]_<COMPONENT>[_ERR]_<OBJECT>_<SHEAR>, plus the three per-object
# metadata columns that carry neither OBJECT nor SHEAR. Anchored so a stray
# legacy token (_PSFo, _Tpsf, NGMIX_ELL_*) fails to match.
# Galaxy is the implicit default object — no token (shapepipe#761), so the
# OBJECT group is OPTIONAL: ``NGMIX_G1_NOSHEAR`` (galaxy, no object token)
# matches alongside ``NGMIX_G1_PSF_ORIG_NOSHEAR`` (explicit PSF object). This
# also makes the regex auto-reject any future ``GAL``-token regression: a
# stray ``_GAL_`` between COMPONENT and SHEAR cannot be absorbed by the
# optional group (``GAL`` is not one of its alternatives) and there is no
# other slot for it to land in.
_COMPONENT = "G1|G2|T|SNR|FLUX|MAG|FLAGS"
_OBJECT = "PSF_ORIG|PSF_RECONV"
_SHEAR = "NOSHEAR|1P|1M|2P|2M"
GRAMMAR_RE = re.compile(
    rf"^NGMIXm?_(?:{_COMPONENT})(?:_ERR)?(?:_(?:{_OBJECT}))?_(?:{_SHEAR})$"
    rf"|^NGMIXm?_(?:MCAL_FLAGS|MCAL_TYPES_FAIL|N_EPOCH|NEIGHBOUR_FLAG)$"
)

# The shipped final-catalogue param files, two levels up from tests/module/.
# Both are consumer contracts updated to the new grammar, so both are checked.
_EXAMPLE = Path(__file__).resolve().parents[2] / "example"
PARAM_PATHS = [
    _EXAMPLE / "cfis" / "final_cat.param",
    _EXAMPLE / "unions_800" / "cat_matched.param",
]

# The one param-file NGMIX token outside the _save_ngmix_data grammar: the
# moments-failure flag, set by a different code path (not the metacal writer).
# Named explicitly so the consistency check stays honest — it is excluded by
# name, not silently dropped, and the test asserts it really is absent from the
# producible set so this carve-out can never mask a real grammar regression.
NON_WRITER_PARAM_TOKENS = {"NGMIX_MOM_FAIL"}


def _param_ngmix_tokens(param_path):
    """Every distinct NGMIX_* token named in a shipped param file."""
    text = param_path.read_text()
    return set(re.findall(r"\bNGMIX[A-Za-z0-9_]*", text))


def _produced_columns(obj_ids):
    """The set of column names the writer emits for these obj_ids."""
    rows = [_base_row(oid) for oid in obj_ids]
    import tempfile
    with tempfile.NamedTemporaryFile(suffix=".fits", delete=False) as fh:
        path = fh.name
    _write_ngmix_cat(path, rows)
    out = _run_save_ngmix(path, obj_ids)
    Path(path).unlink()
    return set(out)


# ---------------------------------------------------------------------------
# Strategies (physically-sensible ranges so generated inputs are valid)
# ---------------------------------------------------------------------------

# Distinct, finite object NUMBERs.
obj_ids_strategy = st.lists(
    st.integers(min_value=1, max_value=10_000),
    min_size=1, max_size=5, unique=True,
)

# A finite, well-separated value: ellipticity-/size-scale magnitudes the writer
# copies verbatim. Bounded away from zero and from each other (see below).
_finite = st.floats(
    min_value=-1.0, max_value=1.0,
    allow_nan=False, allow_infinity=False,
).filter(lambda x: abs(x) > 1e-6)


# ---------------------------------------------------------------------------
# (a) Un-aliasing for ALL inputs
# ---------------------------------------------------------------------------

@settings(deadline=None, max_examples=25)
@given(
    orig_seed=_finite,
    reconv_seed=_finite,
    obj_ids=obj_ids_strategy,
)
def test_psf_orig_reconv_unaliased_for_all_distinct_sources(
    orig_seed, reconv_seed, obj_ids, tmp_path_factory
):
    """Orig/reconv PSF columns trace to their OWN source, never crossed (#749).

    For arbitrary distinct seeds, fill the six orig-PSF source keys with values
    derived from ``orig_seed`` and the six reconv-PSF source keys from
    ``reconv_seed`` — each of the twelve keys a DISTINCT value (seed times a
    per-key multiplier). Then assert every emitted ``*_PSF_ORIG_*`` /
    ``*_PSF_RECONV_*`` column equals the value of the exact source key it must
    read. A family swap (orig<->reconv) OR a within-family swap (G1<->T,
    value<->error) breaks at least one equality, so the property is non-vacuous.
    """
    # Multipliers so each of the six keys per family gets its own value; the
    # two families also never collide because orig_seed != reconv_seed is
    # enforced below and the multiplier set is shared.
    mult = {
        "g1": 1.0, "g2": 1.3, "g1_err": 1.7, "g2_err": 2.1,
        "T": 2.6, "T_err": 3.1,
    }
    # Guard against orig_seed and reconv_seed mapping any key to equal values
    # (would make a cross undetectable for that key). assume() would discard;
    # instead nudge reconv so the families are provably disjoint per key.
    if abs(orig_seed - reconv_seed) < 1e-3:
        reconv_seed = reconv_seed + (0.5 if reconv_seed < 0 else -0.5)

    def _src_val(suffix, family_seed):
        # suffix e.g. "g1", "g1_err", "T"; matched to a multiplier by stripping
        # the "_psf_orig"/"_psf_reconv" tail the caller already removed.
        return family_seed * mult[suffix]

    rows = []
    for oid in obj_ids:
        row = _base_row(oid)
        for col_infix, src_key in ORIG_COLS.items():
            suffix = src_key.replace("_psf_orig", "")
            row[src_key] = _src_val(suffix, orig_seed)
        for col_infix, src_key in RECONV_COLS.items():
            suffix = src_key.replace("_psf_reconv", "")
            row[src_key] = _src_val(suffix, reconv_seed)
        rows.append(row)

    ngmix_path = tmp_path_factory.mktemp("ng") / "ngmix.fits"
    _write_ngmix_cat(ngmix_path, rows)
    out = _run_save_ngmix(ngmix_path, obj_ids)

    n = len(obj_ids)
    for shear in SHEAR_EXTS:
        for col_infix, src_key in ORIG_COLS.items():
            suffix = src_key.replace("_psf_orig", "")
            npt.assert_allclose(
                out[f"NGMIX_{col_infix}_{shear}"],
                [_src_val(suffix, orig_seed)] * n,
                err_msg=f"NGMIX_{col_infix}_{shear} not from {src_key}",
            )
        for col_infix, src_key in RECONV_COLS.items():
            suffix = src_key.replace("_psf_reconv", "")
            npt.assert_allclose(
                out[f"NGMIX_{col_infix}_{shear}"],
                [_src_val(suffix, reconv_seed)] * n,
                err_msg=f"NGMIX_{col_infix}_{shear} not from {src_key}",
            )
        # The two families are genuinely distinct everywhere — the un-aliasing.
        assert not np.allclose(
            out[f"NGMIX_G1_PSF_ORIG_{shear}"],
            out[f"NGMIX_G1_PSF_RECONV_{shear}"],
        )
        assert not np.allclose(
            out[f"NGMIX_T_PSF_ORIG_{shear}"],
            out[f"NGMIX_T_PSF_RECONV_{shear}"],
        )


# ---------------------------------------------------------------------------
# (b) Every emitted column name matches the grammar
# ---------------------------------------------------------------------------

@settings(deadline=None, max_examples=25)
@given(obj_ids=obj_ids_strategy)
def test_emitted_column_names_match_grammar(obj_ids, tmp_path_factory):
    """Every column _save_ngmix_data emits matches the grammar regex.

    The grammar is ``NGMIX[m]_<COMPONENT>[_ERR]_<OBJECT>_<SHEAR>`` plus the
    four metadata columns. The regex is anchored, so a legacy token (``_PSFo``,
    ``_Tpsf``, ``NGMIX_ELL_*``) or a malformed name would fail to match. Run
    over hypothesis-varied object sets to confirm the emitted name set is
    input-independent and always well-formed.
    """
    rows = [_base_row(oid) for oid in obj_ids]
    ngmix_path = tmp_path_factory.mktemp("ng") / "ngmix.fits"
    _write_ngmix_cat(ngmix_path, rows)
    out = _run_save_ngmix(ngmix_path, obj_ids)

    assert out, "writer produced no columns"
    for col in out:
        assert GRAMMAR_RE.match(col), f"off-grammar column emitted: {col!r}"
        # Galaxy is the implicit default object — a regressed GAL token would
        # match neither the optional OBJECT group nor any other slot, but
        # assert it directly so the failure message names the real defect.
        assert "_GAL_" not in col and not col.endswith("_GAL"), col


# ---------------------------------------------------------------------------
# (c) Param-file tokens are all producible by the writer
# ---------------------------------------------------------------------------

@pytest.mark.parametrize(
    "param_path", PARAM_PATHS, ids=lambda p: f"{p.parent.name}/{p.name}"
)
@settings(deadline=None, max_examples=15)
@given(obj_ids=obj_ids_strategy)
def test_param_file_ngmix_tokens_are_producible(param_path, obj_ids):
    """Every NGMIX_* token the param file names is a column the writer produces.

    Each shipped final-catalogue param file (``example/cfis/final_cat.param``
    and ``example/unions_800/cat_matched.param``) is a consumer contract for the
    final catalogue; ``create_final_cat`` keeps only the listed columns, so a
    token it names that the writer cannot emit is a silent, empty column
    downstream. Both files are checked so a future divergence in either (a
    typo'd or stale NGMIX token) cannot escape the consistency check. The one
    known exception — ``NGMIX_MOM_FAIL``, the moments-failure flag set by a
    different path — is excluded BY NAME, and we assert it is genuinely outside
    the producible set so the carve-out cannot hide a real grammar regression.

    Driven over varied object sets to confirm the producible set is stable
    across inputs (the contract is structural, not data-dependent).
    """
    produced = _produced_columns(obj_ids)
    param_tokens = _param_ngmix_tokens(param_path)

    # The carve-out is real: the excluded token is NOT producible. (If a future
    # change made the writer emit it, this fails and we drop the exclusion.)
    for tok in NON_WRITER_PARAM_TOKENS:
        assert tok not in produced, (
            f"{tok} is now producible — remove it from NON_WRITER_PARAM_TOKENS"
        )

    missing = (param_tokens - NON_WRITER_PARAM_TOKENS) - produced
    assert not missing, (
        f"param file names NGMIX columns the writer cannot produce: "
        f"{sorted(missing)}"
    )

    # Non-vacuous: the writer DOES cover real grammar tokens from the param
    # file (so the assertion above is exercised against a non-empty set).
    assert (param_tokens - NON_WRITER_PARAM_TOKENS) & produced


# ---------------------------------------------------------------------------
# (d) The FULL frozen grammar, across all three estimator families
# ---------------------------------------------------------------------------
#
# ``ESTIMATOR_COMPONENT[_ERR]_OBJECT[_metacaltype]`` (shapepipe#761). This
# section encodes the grammar itself as a static property — example column
# names, not a driven writer — so it covers ngmix and HSM together without
# depending on the HSM producer module (covered by its own writer-driving
# test: ``test_psfex_interp.py``). What it pins:
#
# - galaxy is the implicit default object, no token, for the g-type ngmix
#   writer;
# - ellipticity is split into named scalar components, never a packed
#   2-vector: ``G1``/``G2`` for g-type (ngmix, HSM since the return to
#   reduced shear), never a packed ``ELL`` or e-type ``E1``/``E2``;
# - exactly one stored size, ``T`` (no ``SIGMA``/``FWHM``/``r50``);
# - HSM keeps the singular ``FLAG`` token; ngmix keeps plural ``FLAGS``.

FROZEN_GRAMMAR_RE = re.compile(
    r"^(?:"
    # ngmix: galaxy (no object token) or explicit PSF_ORIG/PSF_RECONV.
    r"NGMIXm?_(?:G1|G2|T|SNR|FLUX|MAG|FLAGS)(?:_ERR)?"
    r"(?:_(?:PSF_ORIG|PSF_RECONV))?_(?:NOSHEAR|1P|1M|2P|2M)"
    r"|NGMIXm?_(?:MCAL_FLAGS|MCAL_TYPES_FAIL|N_EPOCH|NEIGHBOUR_FLAG)"
    # HSM: g-type, explicit PSF/STAR object, singular FLAG; the multi-epoch
    # sink in make_cat._save_psf_data appends a bare epoch index.
    r"|HSM_(?:G1|G2|T)_(?:PSF|STAR)(?:_\d+)?"
    r"|HSM_FLAG_(?:PSF|STAR)(?:_\d+)?"
    r")$"
)

# Representative columns the new grammar PRODUCES — one or more per family,
# chosen to exercise every rule above (GAL-drop, G1/G2 vs E1/E2, single T,
# FLAG vs FLAGS).
_GRAMMAR_VALID_EXAMPLES = [
    "NGMIX_G1_NOSHEAR",  # galaxy, no object token
    "NGMIX_T_ERR_PSF_ORIG_1M",
    "NGMIX_FLAGS_2P",
    "NGMIXm_G1_PSF_RECONV_NOSHEAR",
    "NGMIX_MCAL_FLAGS",
    "NGMIX_NEIGHBOUR_FLAG",  # per-object blend flag (shapepipe#776)
    "HSM_G1_PSF",
    "HSM_T_STAR",
    "HSM_FLAG_PSF",
    "HSM_G1_PSF_3",  # make_cat multi-epoch sink
    "HSM_T_PSF_2",
    "HSM_FLAG_STAR_1",
]

# Pre-#761 / off-grammar columns the rename REMOVES — each violates exactly
# one rule, named so a failure here points at the specific regression.
_GRAMMAR_INVALID_EXAMPLES = [
    "NGMIX_G1_GAL_NOSHEAR",  # GAL token regression
    "NGMIX_ELL_NOSHEAR",  # packed ellipticity
    "E1_PSF_HSM",  # pre-rename HSM naming (not HSM_-prefixed)
    "HSM_E1_PSF",  # e-type leftover: HSM stores g-type G1/G2 since #761
    "SIGMA_PSF_HSM",  # raw sigma, not T
    "HSM_SIGMA_PSF",  # stored sigma instead of T
    "HSM_FLAGS_PSF",  # plural — HSM is singular FLAG
    "NGMIX_ELL_PSF_ORIG_NOSHEAR",  # packed ellipticity, not G1/G2
    "SPREAD_MODEL",  # removed entirely, not renamed
]


@pytest.mark.parametrize("col", _GRAMMAR_VALID_EXAMPLES)
def test_frozen_grammar_accepts_new_grammar_examples(col):
    """One representative column per family/rule matches the frozen grammar."""
    assert FROZEN_GRAMMAR_RE.match(col), f"should match grammar: {col!r}"


@pytest.mark.parametrize("col", _GRAMMAR_INVALID_EXAMPLES)
def test_frozen_grammar_rejects_legacy_examples(col):
    """Each pre-#761 / removed column is rejected, not silently absorbed."""
    assert not FROZEN_GRAMMAR_RE.match(col), f"should NOT match grammar: {col!r}"
