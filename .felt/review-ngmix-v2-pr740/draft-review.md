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
| **CI passing** | ⚠️ **Build GREEN, suite NOT run.** Via same-repo mirror PR #741 (pushed `ngmix_v2.0` to origin). The branch's *old* `deploy-image.yml` ran and **passed every step**: builds both images (so the new `uv.lock` + ngmix 2.4.0 resolve and the image builds with all system tools ✅), binary smokes ✅, entry point ✅, `pytest --version` ✅, published to ghcr ✅. But it does **not** run the pytest suite or example pipeline — those need the modern workflow. The modern test-running workflow (`pull_request → develop`) doesn't queue because GitHub evaluates the trigger from the *head branch's* workflow file, which predates the modernization. `ci-release.yml` runs the full suite only on `pull_request → main/master`. **To run the full suite + example pipeline, the branch must have `develop` merged in** (pulls in modern CI + resolves drift). Original #740 got no CI at all because fork PRs don't trigger our Actions without maintainer approval. |
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

## Line-level (staged as draft comments)

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
