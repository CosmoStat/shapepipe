# ngmix v2.0 weight / inverse-variance handling — decision-ready report

**Scope:** PR #741 (branch `ngmix_v2.0`), the v2.0 refactor of noise/weight handling in `ngmix.py`, against the v1 baseline (`origin/develop`), and its relationship to Issue #604 (Weight Handling). All file:line refs are in the worktree (`ngmix_v2.0` head + review cleanups). Load-bearing claims verified against source **and confirmed empirically** (see §2a).

**Headline:** v2.0 introduced two real, test-verifiable regressions in the ngmix weight map (dead noise estimator + lost mask binarization). The minimal v1-restore belongs in #741; the proper inverse-variance fix from #604 is a clean, mostly-config follow-up PR that should land separately.

---

## 1. Precise regression characterization

The entire live weight path is `prepare_ngmix_weights` (`ngmix.py:871–908`):

```python
weight_map = np.copy(weight)        # 894  FSCALE-rescaled Megapipe 0/1 mask
weight_map[flag != 0] = 0.0         # 895
sig_noise = sigma_mad(gal)          # 897  noise scalar over the WHOLE stamp
...
weight_map *= 1 / sig_noise ** 2    # 906
```

v1 (`origin/develop`, `do_ngmix_metacal`), in order:

```python
weight_map = np.copy(weights[n_e])
weight_map[np.where(flags[n_e] != 0)] = 0.0
weight_map[weight_map != 0] = 1     # binarize to exactly 0/1   <-- DROPPED in v2
...
sig_noise = get_noise(...)          # object-free windowed estimator  <-- ORPHANED in v2
...
weight_map *= 1 / sig_noise**2
```

Two coupled regressions, both verified:

**R1 — Noise estimator regressed from object-free to whole-stamp.** v1 called the windowed `get_noise` (masks the object via a model-Gaussian window, re-estimates `sigma_mad` on object-free pixels). v2 uses bare `sigma_mad(gal)` over the full, flux-contaminated stamp. `get_noise` still exists (`ngmix.py:826`, body byte-identical to v1) but **has zero call sites in `src/`** — confirmed dead. v2's `sig_noise` is biased high by the galaxy's own flux, in a **size/flux-dependent** way. Dominant numerical effect for today's inputs.

**R2 — Lost the mask binarization.** v1 collapsed the weight to exactly 0/1 before scaling, so the final weight was a clean uniform `1/sig_noise²` per unmasked pixel. v2 dropped that line, so line 906 multiplies the *raw rescaled mask* by `1/sig_noise²`.

What `weight` contains: traced from `ngmix_runner.py` input #4 → vignetmaker `RUN_2` → single-exposure weight images = Megapipe weight maps, which Gwyn confirmed are **0/1 masks** (#604). The only transform before line 894 is the FSCALE rescale (`rescale_epoch_fluxes`: `weight * 1/Fscale**2`) — identical to v1.

Consequences of R2:
- **For today's 0/1 masks:** each unmasked pixel carries `Fscale⁻² / sig_noise²`. FSCALE is applied to the weight *twice* and to the noise *once* → per-pixel ivar uniformly off by `1/Fscale²` per epoch vs v1. A latent scale bug, varies epoch-to-epoch.
- **Latent hazard (dangerous):** if a *real* inverse-variance map is fed in (THELI, or the #604 BACKGROUND_RMS path), v2 multiplies a genuine `1/var_pix` by another `1/sig_noise²` → **variance counted twice**. v1's binarization made the code robust to this; v2 is not.

### Downstream direction of the effect
`weight_map` → `Observation.weight` = per-pixel inverse-variance in the ngmix χ².
- Inflated `sig_noise` → too-low weight → under-weighted likelihood → **over-estimated** `g_err`/`T_err`/`flux_err`, **under-estimated** `s2n`. Because whole-stamp `sigma_mad` grows with brighter/larger galaxies, this is a **size/flux-dependent miscalibration** — the signature of a multiplicative shear bias `m`, consistent with the prior old-path finding `m ~ -2.8e-2`.
- **Multi-epoch PSF averaging** weights each epoch by `obs.weight.sum()` → now driven by whole-stamp MAD and FSCALE; reported `PSFo` g/T shift accordingly.

---

## 2. Test-verifiability: **YES — clean red→green**

`make_data` (`simulate.py:91`) injects the exact truth inverse-variance: `weights = im*0 + 1/noise**2`, `flags = 0`. That `1/noise²` is the oracle. A unit-level test on `prepare_ngmix_weights` (no ngmix fit needed — fast, deterministic) asserts the returned weight equals the injected inverse-variance:

```python
def test_weight_map_recovers_injected_inverse_variance():
    noise = 1e-3
    gals, psfs, _, weights, flags, jacobs = make_data(
        rng=np.random.RandomState(123), shear=(0.0, 0.0),
        noise=noise, n_epochs=1, img_size=201)
    _, weight_map, _ = prepare_ngmix_weights(
        gals[0], weights[0], flags[0], np.random.RandomState(0))
    truth_ivar = 1.0 / noise ** 2
    recovered = np.median(weight_map[weight_map > 0])
    npt.assert_allclose(recovered, truth_ivar, rtol=0.10)
```

**Coverage gap today:** `test_ngmix.py`'s only weight-path test (`test_metacal_is_reproducible_with_fixed_seed`) asserts `run(42)==run(42)` — a determinism check that passes for *any* deterministic normalization, correct or not. The regression is exercised but never checked against truth.

**Caveats:** `make_data` injects *homoscedastic* noise — the test catches the double-normalization + scalar bias, but to prove the fix respects *spatially-varying* ivar (#604's point), `make_data` needs a per-pixel RMS option (follow-up, not a blocker).

### 2a. EMPIRICAL CONFIRMATION (observed, this session)

Ran the check against the actual code (dev container, galsim+numpy; no ngmix fit needed):

| noise | truth ivar | sigma_mad(gal) | recovered median | ratio (correct = 1.0) |
|---|---|---|---|---|
| 1e-3 | 1e6 | 1.065e-3 | 8.81e11 | **8.81e5** |
| 1e-4 | 1e8 | 1.167e-4 | 7.35e15 | **7.35e7** |

The recovered weight is off by ≈`1/noise²` — exactly the R2 double-count (real ivar × another `1/sig_noise²`). `sigma_mad(gal)` runs ~6% above truth noise (R1 flux contamination, small for a compact galaxy in a 201² stamp). The data matches the by-reading analysis to ~6 orders of magnitude. **The regression is real and observed, not inferred.**

---

## 3. Recommended baseline: SExtractor BACKGROUND_RMS — concrete change plan

The existing machinery is already generic: SExtractor's check-image handler is list-driven and image-agnostic; vignetmaker's ME loop loops over `ME_IMAGE_PATTERN` and reads `get_data(0)` — a single-HDU BACKGROUND_RMS check-image flows through the same path as BACKGROUND with **zero code change**.

| Step | What | File(s) | Cost |
|---|---|---|---|
| 1 | `CHECKIMAGE = BACKGROUND, BACKGROUND_RMS` | `example/cfis/config_exp_*.ini` | **config only** |
| 2 | Add `background_rms` to vignetmaker `ME_IMAGE_DIR`/`ME_IMAGE_PATTERN` | `config_tile_MiViSmVi.ini` | **config only** |
| 3 | Wire one **optional** input (RMS sqlite) through runner + `Vignet` | `ngmix_runner.py`, `ngmix.py` `Vignet`/`Ngmix.__init__` | small code |
| 4 | **Core:** build per-pixel ivar from RMS, gated on mask | `ngmix.py` `prepare_ngmix_weights` et al. | real work |
| 5 | Keep `sigma_mad`/`get_noise` as **fallback** when no RMS map present | `ngmix.py` (branch on `bkg_rms is not None`) | small code |

**Core (Step 4) sketch** — Megapipe weight stays the *mask*; RMS supplies the *variance* (the two roles were conflated before):

```python
weight_map = np.copy(weight)        # Megapipe mask (0 where masked)
weight_map[flag != 0] = 0.0
mask = weight_map != 0
ivar = np.zeros_like(gal)
ivar[mask] = 1.0 / bkg_rms[mask] ** 2          # per-pixel inverse variance
sig_noise = np.median(bkg_rms[mask]) if mask.any() else sigma_mad(gal)
return gal_masked, ivar, noise_img
```

RMS rescaling: BACKGROUND_RMS is in image counts, scales by `Fscale`; `1/rms²` scales by `1/Fscale²` — identical to the existing weight rescale, slots into `rescale_epoch_fluxes`.

**Effort: S–M, ~4–8 h.** Config + plumbing ~1–2 h; core rewrite + fallback ~2–3 h; example-pipeline before/after on shear/s2n/flags ~2–3 h.

**Honest read on aguinot's "almost free":** true for the plumbing (steps 1–2 are pure config), undersells the core. Delivering the data is cheap; *using it correctly* (changing what the weight means, getting Fscale + mask roles right, validating no fit regression) is the earned part.

**Risks:** RMS units/rescaling; zeros/NaNs in RMS at masked/off-image pixels → `1/rms²` blowup (gate on mask, guard `rms<=0`); double-counting the mask; positional input-list shift breaking configs (mitigated by optional-input design).

---

## 4. Scoping: **separate PRs.**

- **In #741 (now):** the minimal v1-restore — reinstate `weight_map[weight_map != 0] = 1` and call `get_noise(...)` instead of `sigma_mad(gal)` — plus the red→green unit test. Two-line change fixing a regression *this PR introduced*; restores known-good behavior; restoring the binarization makes the fallback robust to real ivar maps again. Does not depend on #604.
- **#604 as its own PR (next):** the BACKGROUND_RMS baseline changes what the weight map *means*, touches new I/O, needs its own before/after validation. A distinct review surface; bundling into #741 would mix "fix the regression I just introduced" with "implement a 3-year-old feature request."

The red→green test validates both the #741 minimal fix and the #604 baseline; it travels with #741 now and proves out #604 later.

---

## 5. Open questions for Martin (science call)

1. **Baseline vs THELI.** Confirm in-house SExtractor-BACKGROUND_RMS is the baseline, with THELI weight images slotting in later as another `ME_IMAGE_PATTERN` source (not blocking).
2. **Optional vs required RMS input.** Recommend optional with `sigma_mad`/`get_noise` fallback. Confirm.
3. **#741 minimal fix:** restore the windowed `get_noise` (reintroduces the gal-guess machinery the PR may have deliberately stripped), or accept bare `sigma_mad(gal)` until #604 lands? R1 is the dominant effect, so this matters even short-term.
4. **Does the RMS map feed anything beyond the weight + noise-fill scale** (e.g. background subtraction)? Recommend weight + noise-fill only.
5. **Scope of config edits** for the #604 PR: cfis only, or cfis + cfis_simu + canfar batch?
6. **Validation bar to merge #604:** smoke before/after on the example tile, or a full sim m/c calibration before merge?
