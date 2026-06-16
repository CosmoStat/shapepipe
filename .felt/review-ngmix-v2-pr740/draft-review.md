# Draft review — PR #740 "Ngmix v2.0" (staging, NOT posted)

Staged for the call with Cail. Posting/pushing happens *with* him, in his voice.
Authorship convention: "— Claude on behalf of Cail" unless he says otherwise.

---

## Reviewer checklist (Martin's template)

| Box | Verdict |
|---|---|
| Targets `develop` | ✅ head `ngmix_v2.0` → base `develop` |
| Assigned to developer | ⬜ (Cail to check on GitHub) |
| Appropriate labels | ☑️ `enhancement` only — fine; could add a `dependencies` label given the ngmix bump |
| Projects/milestones | ⬜ (Cail to check) |
| Clear description | ⚠️ Summary is terse (4 bullets). For a +3894/−3410, 67-file PR that bumps a core dep and rewrites the shape-measurement module, a fuller description would help — esp. the centroid-bias result. |
| "closes #" links | ❌ No issue links. #725 (Axel's centroid-shift PR) overlaps this work — should this PR close/supersede it? Worth an explicit cross-ref. |
| Code/doc style | ⚠️ Mostly good. Pockets of debug cruft (see line-level). |
| Docs updated | ❌ No `docs/` changes. The ngmix pin note in `CLAUDE.md` ("don't modernize this line") is now stale — this PR *is* the modernization. Needs updating on/with merge. |
| **CI passing** | ✅ **GREEN** (full suite). Got here by pushing a same-repo mirror (PR #741, branch `ngmix_v2.0` on origin) and merging `develop` into it (clean except `uv.lock`, regenerated via `uv lock`) to pull in the modernized CI. `build-test-publish` ran every step green: images build (new `uv.lock` + ngmix 2.4.0 resolve), binary smokes, entry point, **pytest suite** (incl. develop's `test_ngmix.py`), example pipeline; publish correctly skipped on the PR. Original fork PR #740 got no CI (Actions approval gate). | The modern test-running workflow (`pull_request → develop`) doesn't queue because GitHub evaluates the trigger from the *head branch's* workflow file, which predates the modernization. `ci-release.yml` runs the full suite only on `pull_request → main/master`. **To run the full suite + example pipeline, the branch must have `develop` merged in** (pulls in modern CI + resolves drift). Original #740 got no CI at all because fork PRs don't trigger our Actions without maintainer approval. |
| API docs built | ➖ n/a-ish (no public API doc surface changed materially) |
| All files checked + comments | ✅ done in this review |

---

## Science-quality read (the part the checklist misses)

### 1. ngmix 2.4.0 API correctness — ✅ verified empirically
Installed ngmix **2.4.0** (from `esheldon/ngmix@v2.4.0`, the PR's own source) over the
container stack and ran the **actual module code**. All 13 ngmix symbols the module
uses resolve in 2.4.0 (`fitting.Fitter`, `guessers.TPSFFluxAndPriorGuesser`/`TFluxGuesser`,
`runners.PSFRunner`/`Runner`, `metacal.MetacalBootstrapper`, `joint_prior.PriorSimpleSep`,
`priors.*`, `observation.*`, `Jacobian`). The metacal/bootstrap path runs end-to-end.

### 2. Shear recovery — ✅ unbiased
Ran `do_ngmix_metacal` on 200 galsim galaxies (sub-pixel shifted, the centroid regime),
recovered shear via the metacal response:
- **m = +2×10⁻⁴ ± 3×10⁻⁴** (multiplicative bias consistent with zero)
- **c₂ = −1×10⁻⁵** (additive consistent with zero), R₁₁≈R₂₂≈0.924, 200/200 fits OK.
- Author's own `centroid_bias_v2.py` (noiseless, 60 trials, ±-cancellation): **m = +3.7×10⁻⁴ ± 0.5×10⁻⁴ (1σ)**, c₂ = (−0.3 ± 0.15)×10⁻⁵, R₁₁≈0.9234.
- **Read:** few×10⁻⁴ residual m — small, likely within UNIONS requirements, but the noiseless limit resolves a marginal *positive* residual (≈7σ given tiny noiseless bars). Not the centroid fix (fix-off gives the same) → metacal/fitter higher-order. Worth stating the m-requirement we're holding this to.

### 3. The centroid fix — present, correct in form, but **necessity not demonstrated in idealized sim**
The fix (`make_ngmix_observation`, `ngmix.py:987-1004`) re-centers the galaxy Jacobian on
the HSM `FindAdaptiveMom` centroid so the centroid prior (centered at the Jacobian origin)
doesn't bias the fit. Sound idea. **But** running the same ensemble with the recenter
disabled gives the *same* unbiased m (both ~3×10⁻⁴) — i.e. in this small-shift (±½ px),
high-S/N regime the fix is harmless but not measurably necessary. **Question for the call:**
in what regime does the bias the fix targets actually appear? The author's centroid scripts
(`centroid_bias.py`, configs bug/fix/v2) presumably demonstrate it — is the before/after
result captured anywhere we can point to? This is the single most important thing to pin down.

### 4. `r50` vs `T` labeling — ⚠️ likely mislabel, worth clarifying
In `compile_results` (`ngmix.py:415-417`): `output_dict[name]['r50'] = results[...]['pars'][4]`.
For an ngmix `gauss` model, `pars[4]` is **T** (≈2σ², arcsec²), *not* a half-light radius —
and `output_dict[name]['T']` stores the same value. And `r50_psfo = sqrt(T_psf/2)` is the
Gaussian **σ**, not r50 (true Gaussian r50 = 1.1774σ). So columns named `r50*` actually carry
T / σ. Downstream (sp_validation) needs to know what these columns mean. Is "r50" intentional
shorthand, or should these be true half-light radii? (Constitution flagged "r50 as the size
guess changed from T" — but the *guesser* still uses `T=0.25`; the change is in output column
naming, not the guess. Worth confirming Martin's intent.)

---

## Code-review pass (multi-agent, high-effort) — findings

Ran a 6-angle code-review over the ngmix delta. CI is green, so none of these break the
test suite — but the suite is shallow on these paths. Ranked.

### Correctness — confirmed
1. **`get_guess` ImportError — shipped script dead on arrival.** `scripts/validation/centroid/centroid_bias.py:37` does `from ...ngmix import get_guess, get_noise`, but this PR removed `get_guess` from the module. The script ImportErrors on load. (Only `centroid_bias.py`; `test_centroid_shift.py` imports just `get_noise`.) Fix the import, restore `get_guess`, or drop the v1 script.
2. **Reproducibility break — unseeded RNG.** `ngmix.py:948-949` `prepare_ngmix_weights` uses global `np.random.randn` for the noise image + masked-pixel fill, while the module deliberately switched to a seeded `self._rng` (`:266`). Same tile + same seed → different masked-pixel values and metacal noise → non-reproducible shears. Thread `rng` into `prepare_ngmix_weights`/`make_ngmix_observation`.
3. **`*_psfo` columns now carry the metacal-RECONVOLVED PSF, not the original PSFEx/MCCD PSF.** `average_multiepoch_psf` reads `obsdict['noshear'][n_e].psf.meta['result']` — the reconvolved metacal PSF — and feeds `g_PSFo`/`T_PSFo`/`r50_PSFo`. But `compile_results` documents these as "the original image psf from psfex or mccd," and PSF-leakage / ρ-stat / null tests consume them as the input PSF. The old code fit the raw pre-metacal PSF separately. Confirm intent; the comment is now wrong either way.
4. **`CHECK_EXISTING_DIR` resume silently dropped.** New `process()` never reads `self._check_existing_dir` or calls `get_last_id` (both now dead). A configured resume reprocesses from scratch and overwrites — a documented feature silently no-ops.

### Correctness — plausible, needs Martin/Lucy intent (raise as questions)
5. **Weight-map normalization changed.** Old `prepare_ngmix_weights` binarized the mask (`weight[weight!=0]=1`) then applied a flat `1/sig_noise²`. New keeps the per-pixel (already FSCALE-rescaled) weight AND multiplies by `1/sigma_mad(gal)²`, with `sigma_mad` taken over the *object-containing* stamp (flux-contaminated, vs the old windowed object-free `get_noise`). Looks like double inverse-variance + a biased noise estimate → shifts the fit's χ² weighting, errors, s2n, and PSF-averaging weights vs the validated v1. Intentional?
6. **Per-type PSF columns collapsed.** `Tpsf`/`g1_psf`/`g2_psf` are now identical across all 5 metacal HDUs (all from `T_PSFo`/`g_PSFo`); old code stored per-type values. May be fine if the metacal target PSF is type-independent — but if any response/selection step differences them it now gets exactly zero.
7. **Galaxy zero-pixel guard narrowed (`any`→`all`).** `prepare_postage_stamps:737` now skips an epoch only if `np.all(gal_vign==0)`; the old code skipped if *any* pixel was 0. Partial CCD-edge stamps (a band of exact-zero pixels) now enter the fit, relying on flags/weights to catch them. Intentional loosening or a regression?
8. **PSF fit uses the galaxy joint prior + galaxy `FLUX_AUTO` as the PSF flux guess** (`do_ngmix_metacal`). Matches the upstream ngmix example's prior reuse, but the galaxy-flux guess for the PSF is odd; possible PSF-fit convergence/quality impact. Lower confidence.

### Altitude / structure
9. **Runner input-contract mismatch (4 runners).** `sextractor_runner`, `read_ext_sexcat_runner`, `psfex_interp_runner`, `vignetmaker_runner` moved the WCS/exp-numbers paths from named config options (`LOG_WCS`/`ME_LOG_WCS`/`ME_DOT_PSF_DIR`) to positional `input_file_list[i]` — but, unlike `ngmix_runner`, their `@module_runner` decorators were **not** updated. The real input contract now lives entirely in each config's `FILE_PATTERN`/`FILE_EXT` ordering; the shipped v2.0 example configs are correct, but a hand-written/v1-style config silently mis-maps inputs or IndexErrors with no schema to catch it. `ngmix_runner` is the model — extend the other four decorators to match.
10. **`r50`/`T` mislabel (units).** `compile_results:415-417` stores `pars[4]` (= `T`, arcsec²) under `r50`/`r50_err`, and `r50_PSFo = sqrt(T_psf/2)` is the Gaussian σ, not the half-light radius (`r50 = 1.1774σ`). sp_validation consumes these size columns. Intentional shorthand or a mislabel?

### Cleanups (bundle into one comment)
- `save_results`: `ngmix.py:547` non-`f` string → error prints literal `{ext_name}`; leftover `print(...)` debug at `:564-565`; the `np.append`-per-batch FITS append is O(n²) over a tile.
- Dead code: unreachable `print` after `return` in `MegaCamFlip` (`:293`); unused `sextractor_e1e2` (`:1131`, with a copy-pasted wrong docstring); `# I still don't know how to handle this` (`:63`); commented `#self._gal_vignet_path` block (`:240-244`); `# MKDEBUG` markers; redundant `self._f_wcs_path`.
- `51 * 51` hardcoded masked-fraction denominator (`:766`) → use `v_flag_tmp.size`; near-term goal is configurable stamp sizes, so this will silently go wrong.
- Duplication: `scripts/python/fitting.py` and `scripts/jupyter/test_centroid_shift.py` each redefine their own `make_data`/`get_prior`, duplicating `testing/simulate.py` + the module `get_prior` added in this same PR (defeats simulate.py's "stable across branches" purpose). `fitting.py` is a stray upstream-ngmix example (`fitting_bd_empsf.py` header, `espy`/`images` refs) — likely an unintended add; delete or make it import the canonical helpers.
- Minor efficiency: doubled `galsim.Image(gal, scale=1.0)` in `make_ngmix_observation`; repeated uncached `SqliteDict[str(obj_id)]` deserialization per object/epoch; `get_exp_output_dirs` re-globs identical dirs when a runner repeats in `ME_IMAGE_EXP_RUNNERS`.

---

## Line-level (earlier staged notes — superseded by the code-review pass above for ngmix.py)

- **`pyproject.toml`** — ngmix now `>=2.4` with `[tool.uv.sources] ngmix = git esheldon/ngmix tag v2.4.0`. Why pin to a git tag rather than the PyPI release `ngmix==2.4.0`? PyPI would be more reproducible and drop the git dependency. (PyPI was unreachable from here to confirm 2.4.0 is published — Cail can check.)
- **`ngmix.py:63`** — stray `# I still don't know how to handle this` above `Tile_cat`.
- **`ngmix.py:254`,`:528`-ish** — `# MKDEBUG` markers and `print(...)` debug statements (e.g. `:293` unreachable `print('rotating megapipe image')` after a `return`; `:564-565` `print("output_dict"...)` in an error branch). Clean up before merge.
- **`ngmix.py:766`** — masked-fraction cut hardcodes `51 * 51`; if vignet size ever changes this silently breaks. Use `flag_vign.size`.
- **`scripts/python/fitting.py`** — looks like a lightly-edited copy of esheldon's `fitting_bd_empsf.py` example; docstring still references `fracdev` / `fitting_bd_empsf.py` and a bulge+disk fit that no longer applies. Either trim to match or drop — is it meant to ship?
- **`scripts/jupyter/test_centroid_shift.py`** — jupyter-exported (cell markers `# In[14]:`, commented-out code, inline-duplicated `make_data`). Now that `shapepipe.testing.simulate.make_data` exists, this could import it instead of carrying its own copy. Ship as-is or tidy?
- **`src/shapepipe/testing/simulate.py`** — 👍 clean, dependency-light (galsim+numpy only), gives exactly the cheap validation path. Nice addition.
- **Supporting runners** (`psfex_interp_runner`, `vignetmaker_runner`, `ngmix_runner`) — the `$SP_EXP`-tree resolver (`get_exp_output_dirs`) with config-gated fallback to the v1 run-log path is clean and backward-compatible. `f_wcs_path` now threaded via `input_file_list` consistently. No concerns.

---

## Suggested PR-level summary comment (draft)
> Reviewed the diff and exercised it on candide against **ngmix 2.4.0** (installed from the
> PR's own source). The 2.x API is used correctly and the metacal path runs end-to-end; on
> simulated galaxies the v2.0 fitter recovers shear with **m ≈ 2×10⁻⁴** (consistent with zero).
> Two things before I can tick the checklist: (1) **CI hasn't run** — it's a fork PR, needs an
> "Approve and run" or a push to origin; (2) I couldn't reproduce a centroid bias in the
> idealized sim with the fix *off*, so I'd like to understand the regime where the centroid fix
> matters — is the before/after result from `centroid_bias.py` captured somewhere? Plus a few
> cleanups (debug prints, the `r50`/`T` column naming, stale `CLAUDE.md` ngmix note). Details inline.
> — Claude on behalf of Cail
