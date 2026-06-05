# ngmix v2.0 size-column naming: ground truth, official convention, recommendation

**Scope:** the five new `r50*` size columns added in ngmix v2.0 (PR #741), what they actually contain, what the official UNIONS papers report, and how to make the producer (shapepipe) → shared (cs_util) → consumer (sp_validation) stack honest and consistent.

**Verification:** code claims verified in the worktree; conversion arithmetic verified numerically; paper convention from arXiv:2605.13549 and arXiv:2204.04798; downstream consumption from the sp_validation source tree.

**Headline:** the v2.0 `r50` columns are mislabeled — galaxy `r50` == the area `T`; PSF `r50psf` == `σ` (not `1.1774σ`). Cleanest fix: make the ngmix writer emit honest, correctly-converted columns at the source (true half-light radius `r50 = 1.1774·√(T/2)` for both galaxy and PSF, alongside the existing `T`), rather than a rename or a downstream-only transform. A *separate* load-bearing bug lives in sp_validation's `T_to_fwhm`, which a clean r50/σ surface lets us retire.

---

## 1. Code ground truth — what each size column actually holds

In ngmix `gauss`, `pars = [cen1, cen2, g1, g2, T, flux]` and **`pars[4] = T = 2·σ²` is an AREA** (arcsec²). The half-light radius is `r50 = 1.1774·σ = 1.1774·√(T/2)`. The `σ = √(T/2)` mapping is confirmed twice in-file (`galsim.Gaussian(sigma=np.sqrt(guess[4]/2))`; `r50_psfo = np.sqrt(max(psf_res['T_psf'],0)/2)`).

| Column (FITS) | Source | Value in σ | True r50 = 1.1774σ? | Status |
|---|---|---|---|---|
| `T`, `T_err` | `results[…]["T"]` | `2σ²` (area) | — | Correctly named area |
| `r50`, `r50_err` | `pars[4]`, `pars_err[4]` | **`2σ²` (area)** | **No** | **MISLABELED: == `T`/`T_err` bit-for-bit. Neither √ nor ×1.1774 applied.** |
| `Tpsf` | `T_PSFo` | `2σ_psf²` (area) | — | Correctly named area |
| `T_psfo_ngmix`, `T_err_psfo_ngmix` | `T_PSFo`, `T_err_PSFo` | `2σ_psf²` | — | Correct; `T_psfo_ngmix` **duplicates** `Tpsf` |
| `r50psf` | `r50_PSFo` = `√(T_psf/2)` | **`σ_psf` (length)** | **No — it's σ, off by 1.1774×** | Genuine length, missing the factor |
| `r50_psfo_ngmix` | `r50_PSFo` | `σ_psf` | No — it's σ | **Duplicates `r50psf`** |
| `r50_err_psfo_ngmix` | `T_psf_err/(2·r50_psfo)` | `σ_psf_err` | No — error on σ | Correct error-prop of σ; NaN when σ=0 |

**Headline hazard:** the same name root `r50` means an **area** on the galaxy side and a **length** on the PSF side. A user ratioing galaxy `r50` against `r50psf` divides an area by a length. **Zero columns in the file are a true half-light radius** — every "r50" is either `T` (galaxy) or `σ` (PSF).

**Arithmetic check (verified):** `√(2·ln2) = 1.17741`; `2.355/1.17741 = 2.000` (why the downstream `T_to_fwhm` is dimensionally wrong — §4).

**v1 provenance:** `origin/develop` wrote sizes **only as `T`**. **All five `r50*` columns are NEW in v2.0**, and the galaxy two were wired to the area `pars[4]` — introduced as mislabeled duplicates.

---

## 2. Official convention + downstream consumption

**UNIONS-3500 WL I — A Galaxy Shape Catalogue (arXiv:2605.13549), the current shape-catalogue paper:**
- Reports the **half-light radius `r_h` (= ngmix `r50`)** as the **PRIMARY** size for **both galaxy and PSF**.
- Size cut is the dimensionless **`0.707 ≤ r_h/r_h,psf ≤ 3`**, applied *inside Metacalibration* (enters the selection response).
- Resolution factor `R = r_PSF/(r_PSF + r_h)` — written in `r`, not `T`.
- **`T` is derived only**, via Eq. 1 `T = (r_h²/ln2)·(1+g1²+g2²)/(1−g1²−g2²)`, "to relate r_h to the DES size definition." Round-object case `T = r_h²/ln2 = 2σ²` — identical to standard ngmix `T = 2σ²`.

**Guinot et al. 2022 (arXiv:2204.04798):** mixed — galaxy cut in `T` (`T_gal/T_PSF > 0.5`), PSF size `T = 2σ²`, but ngmix fit param + a prior in `r50` (arcsec).

**sp_validation:** treats galaxy/PSF size only as the dimensionless **ratio** `T/Tpsf` (size cut, DES-weight & leakage binning) — robust regardless of absolute meaning — plus one place where PSF `T` is converted to FWHM for the leakage fit (§4). **It reads NONE of the five `r50*` columns** (0 hits in `src/`).

**Takeaway:** the official paper wants `r_h` (= r50) as the headline size — exactly the quantity the v2.0 `r50*` columns *claim* to provide but don't. Standardize on **r50 = half-light radius as primary, T = 2σ² as the derived DES area.**

---

## 3. Recommendation — TRANSFORM at the source, thin cs_util surface

Three options: a **pure RENAME** (`r50`→`T`/`σ`) is honest but discards the chance to ship the paper's headline quantity; a **cs_util-only TRANSFORM** (fix on read) leaves a foot-gun baked into every FITS file; **TRANSFORM at the source** (recommended) makes the writer emit honest, paper-consistent columns with cs_util holding the convention once.

### 3a. ngmix.py — every size column honest and on the same scale
Keep all `T*` columns unchanged (correct areas). Change the radius columns to true half-light radii:
- **Galaxy:** `r50 = 1.1774·√(T/2)` (= `√(ln2·T)`); error `r50_err = (1.1774/(2√2))·T_err/√T`. *(Currently append the area `pars[4]`/`pars_err[4]`.)*
- **PSF:** multiply existing σ values by `1.17741` so `r50psf`/`r50_psfo_ngmix`/`r50_err_psfo_ngmix` are true half-light radii, not σ.

After this, **every "r50" column is a length meaning the same thing (1.1774σ)** on both sides, matching UNIONS-3500 I Eq. 1 for round objects.

### 3b. De-duplicate
`T_psfo_ngmix` duplicates `Tpsf`; `r50_psfo_ngmix` duplicates `r50psf`. Retire the redundant names (or document as aliases). Cleanest end state — one name per quantity: galaxy `{T, T_err, r50, r50_err}`, PSF `{Tpsf, Tpsf_err, r50psf, r50psf_err}`. Also fix the asymmetry: a `T_err_psfo_ngmix` exists but no `Tpsf_err` twin.

### 3c. cs_util — single source of truth for the size web
Expose thin, correctly-named, tested converters (home: `cs_util`, since `sp_validation/galaxy.py` already does `from cs_util import cfis`, and cs_util is the shared layer):

```
T_to_r50(T)      = 1.1774 * sqrt(T / 2)      # = sqrt(ln2 * T)
r50_to_T(r50)    = 2 * (r50 / 1.1774)**2
sigma_to_fwhm(s) = 2.355 * s                 # already correct in galaxy.py
T_to_fwhm(T)     = 2.355 * sqrt(T / 2)       # CORRECTED (see §4)
```

---

## 4. Migration impact on sp_validation

**No change for the size-ratio paths.** Galaxy `T` and PSF `Tpsf` enter only as the dimensionless ratio `T/Tpsf` (size cut, `size_ratio` binning) — robust by construction. None of the five `r50*` columns is read in sp_validation `src/` (0 hits), so the mislabel harms no current computation.

**One genuine, load-bearing fix — the PSF-leakage size regressor.** `fwhm_PSF = T_to_fwhm(T_PSFo)` (built in `notebooks/extract_info.py`) is the third regressor in the scale-dependent PSF-leakage fit (`run_object.py:510`, `size_PSF_col` at `calibration.py:540`). The current `galaxy.py:T_to_fwhm` is `T/1.17741*2.355` and its own `MKDEBUG` comment doubts it ("Why is FWHM not quadratic in T???"). It is **wrong**: `T = 2σ²` is an area, so `FWHM = 2.355·σ = 2.355·√(T/2)` — a √ is required, not a linear factor (`2.355/1.17741 = 2.000` shows it treats T as if already σ). Result: a monotonic-but-nonlinear distortion of the PSF-size axis, biasing any `α(PSF-size)` trend and per-bin leakage coefficients. Spatially-constant leakage unaffected.

**Migration:** point the `fwhm_PSF` builder at the new corrected PSF radius column via `cs_util.sigma_to_fwhm` (or import the fixed `cs_util.T_to_fwhm`), retire the local `galaxy.py:T_to_fwhm`. This is the change that actually moves a science number; the column rename/transform is hygiene + paper-consistency that also removes the foot-gun.

**Sequencing:** the catalogue-writer change (§3a/3b) and the sp_validation regressor fix (§4) are independent — sp_validation doesn't read the `r50*` columns today, so fixing ngmix won't break it. cs_util converters should land before sp_validation switches to importing them.

---

## 5. Open questions for Lucy / Martin

1. **Primary on-disk convention:** standardize on **r50 as primary** (matching UNIONS-3500 I) while keeping `T`? Or keep `T` as the working column with r50 derived?
2. **Drop vs keep duplicates** (`T_psfo_ngmix`≡`Tpsf`, `r50_psfo_ngmix`≡`r50psf`): anything outside sp_validation `src/` (notebooks, older catalogues) bound to the `_psfo_ngmix` names? If so, alias rather than delete.
3. **cs_util as the home** for the converter web — agreed? Naming preference (`T_to_r50` vs `T_to_halflight`)?
4. **PSF-leakage refit:** were current/in-flight scale-dependent leakage results produced with the buggy `T_to_fwhm`, and should they be regenerated once fixed? How much does `α(size)` move?
5. **`r50_err` for the galaxy:** is the error-propagation form acceptable, or prefer the full pars covariance?
6. **Catalogue versioning:** does changing on-disk `r50*` semantics warrant a format version bump / changelog so downstream users know v2.0 pre- vs post-fix differ?
