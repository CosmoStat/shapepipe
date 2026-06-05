---
id: 01KTCQPE4M6S1E4V9WMPCAFCKT
name: 'ngmix size columns: honest r50 + cs_util converter web'
tags:
  - constitution
  - shapepipe
  - ngmix
created-at: 2026-06-05T22:30:50.004535516+02:00
outcome: Draft — not yet dispatched.
shuttle:
  enabled: false
  kind: oneshot
  host: candide
  agent: claude-opus
---

# ngmix size columns: honest r50 + cs_util converter web

Spun out of the [[review-ngmix-v2-pr740]] review. The full investigation — code ground-truth table,
the UNIONS-paper convention, the sp_validation consumption map, and the bonus downstream bug — lives in
**`.felt/review-ngmix-v2-pr740/size-report.md`**. Read it first; it is the detailed spec.

The v2.0 ngmix rewrite added five `r50*` columns, **none of which is a half-light radius**: galaxy
`r50` = `pars[4]` = `T` = 2σ² (an *area*, a bit-for-bit copy of the `T` column); PSF `r50psf` =
`√(T/2)` = σ (a length, but off by the `1.1774` factor). So the same name root means an area on the
galaxy side and a length on the PSF side. The official catalogue paper (**UNIONS-3500 WL I**,
arXiv:2605.13549) reports the half-light radius `r_h` (= r50) as the **primary** size for both galaxy
and PSF, with `T` derived. Cail approved fixing this at the source.

## Desired State

**1. ngmix writer emits honest, self-consistent size columns** (`ngmix.py`, on `ngmix_v2.0`):
- Galaxy `r50 = 1.1774·√(T/2)` (= `√(ln2·T)`) with propagated error `r50_err = (1.1774/(2√2))·T_err/√T`
  — replacing the current raw `pars[4]`/`pars_err[4]`.
- PSF `r50psf`/`r50_psfo_ngmix`/`r50_err_psfo_ngmix` = existing σ values × `1.17741` (true half-light
  radii, not σ).
- Keep all `T*` columns as the (correctly-named) DES areas. After this, every "r50" column is a length
  meaning the same thing on both sides, matching the paper.
- **De-duplicate:** `T_psfo_ngmix` ≡ `Tpsf` and `r50_psfo_ngmix` ≡ `r50psf`. Default = retire the
  redundant names; fall back to documenting them as aliases if anything outside `sp_validation/src`
  (notebooks, older catalogues) depends on the `_psfo_ngmix` names. Fix the `T_err_psfo_ngmix`-with-no-
  `Tpsf_err` asymmetry while here.

**2. `cs_util` holds the conversion web** — the single source of truth, since `sp_validation/galaxy.py`
already does `from cs_util import cfis` and cs_util is the shared producer/consumer layer:
- `T_to_r50(T) = 1.1774·√(T/2)`, `r50_to_T(r50) = 2·(r50/1.1774)²`, `sigma_to_fwhm(s) = 2.355·s`,
  and the **corrected** `T_to_fwhm(T) = 2.355·√(T/2)`. With tests.

**3. Fix the downstream PSF-leakage bug in sp_validation.** `galaxy.py:T_to_fwhm` is `T/1.17741*2.355`
— dimensionally wrong (treats the area `T` as σ; `2.355/1.17741 = 2.000`). It builds the `fwhm_PSF`
regressor for the scale-dependent PSF-leakage fit (`run_object.py:510`), so it distorts the PSF-size
axis and biases `α(size)`. Switch it to `cs_util.sigma_to_fwhm` on the corrected PSF size (or import
the fixed `cs_util.T_to_fwhm`) and retire the local function. **This is the change that actually moves
a science number** — flag it prominently for Cail (it touches the pure_eb / leakage analysis).

**Scoping:** the ngmix-writer change (1) and the sp_validation fix (3) are independent — sp_validation
reads none of the `r50*` columns today, so fixing ngmix won't break it. Land cs_util (2) before
sp_validation (3) switches to importing it. **Merge/PR-creation stay Cail's gesture**; surface the open
calls (r50-as-primary convention, drop-vs-alias, catalogue version bump, whether in-flight leakage
results need regenerating) in this fiber's report for Cail's eventual reply — **do not post to Martin.**

## Context

- **Detailed spec:** `.felt/review-ngmix-v2-pr740/size-report.md` (ground-truth table, paper citations,
  consumption map, the 6 open questions).
- **Repos (all on candide):** shapepipe `/automnt/n17data/cdaley/unions/shapepipe` (branch `ngmix_v2.0`);
  sp_validation `/automnt/n17data/cdaley/unions/pure_eb/code/sp_validation/src/sp_validation/`
  (`galaxy.py`, `calibration.py`, `run_object.py`, `notebooks/extract_info.py`); cs_util (shared —
  locate the checkout; `sp_validation/galaxy.py` imports it). Cross-repo, so likely cross-PR.
- **Coordinate with [[ngmix-weights-ivar]]:** both edit `ngmix.py` but different functions (size-column
  writer vs `prepare_ngmix_weights`) — low conflict; each on its own branch off `ngmix_v2.0` (after the
  review-cleanup commit `bd60dc8e` lands).
- **Convention:** UNIONS-3500 WL I (arXiv:2605.13549) → `r_h`/r50 primary, `T = 2σ²` derived (Eq. 1).
  Arithmetic verified: `√(2·ln2) = 1.17741`.
- **Agent:** Cail's call — the cs_util/sp_validation API + naming is taste-y (leans Claude); the ngmix
  column edits are mechanical.
