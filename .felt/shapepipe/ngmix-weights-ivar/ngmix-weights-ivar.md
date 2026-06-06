---
id: 01KTCQPE3JGEYN7NQS8HW1AT6B
name: 'ngmix weight map: fix v2.0 regressions + inverse-variance (#604)'
tags:
  - constitution
  - shapepipe
  - ngmix
created-at: 2026-06-05T22:30:49.970813955+02:00
outcome: |-
  R2 minimal fix landed in `/tmp/pr740-wt` commit `953a52a3`: `prepare_ngmix_weights` again binarizes the mask before scalar noise normalization, and a unit test now goes red on double-normalization (`8.809e11` vs `1e6`) and green after the fix (`8.809e5`, within 15%; residual is R1 whole-stamp MAD contamination).
shuttle:
  enabled: true
  kind: oneshot
  host: candide
  agent: codex
status: active
---

# ngmix weight map: fix v2.0 regressions + inverse-variance (#604)

Spun out of the [[review-ngmix-v2-pr740]] review. The full investigation, with file:line
evidence, the empirical confirmation, the concrete change plan, effort, and risks, lives in
**`.felt/review-ngmix-v2-pr740/weights-report.md`** — read it first; it is the detailed spec.

The v2.0 ngmix rewrite (PR #741, branch `ngmix_v2.0`) introduced two coupled regressions in
`prepare_ngmix_weights` (`ngmix.py:871`), **empirically confirmed** this session (fed `make_data`'s
truth inverse-variance `1/noise²` → recovered `8.8e11`, off by ≈`1/noise²`):

- **R1** — noise estimate regressed from v1's object-free windowed `get_noise` (still present at
  `:826`, now **dead**) to flux-contaminated whole-stamp `sigma_mad(gal)`. Size/flux-dependent bias
  on the likelihood weighting → biased `g_err`/`T_err`/`s2n`; the fingerprint of a multiplicative
  shear bias (cf. old-path `m≈-2.8e-2`). Dominant today.
- **R2** — lost the v1 binarization `weight_map[weight_map != 0] = 1`, so the raw weight is multiplied
  by `1/sig_noise²`: a per-epoch `1/Fscale²` error today, and a latent double-count the moment a real
  inverse-variance map is fed in.

Neither v1 nor v2 is *correct* in absolute terms — both approximate per-pixel inverse-variance with a
per-stamp scalar on a 0/1 mask. #604 ("Weight Handling", open since 2023) is the real fix: feed ngmix
genuine inverse-variance maps.

## Desired State

**Two PRs, in this order:**

1. **Minimal regression fix on `ngmix_v2.0` (#741).** Restore the v1 binarization (R2) and add a
   **red→green unit test** on `prepare_ngmix_weights` using `make_data`'s truth ivar (the test sketch
   is in the report; it needs no ngmix fit, so it's fast and deterministic). This is a clean, small
   change that removes the double-count hazard.
   - **Open call (default = defer R1):** restoring the windowed `get_noise` for R1 reintroduces the
     gal-guess machinery this PR may have deliberately stripped. Default to fixing R2 + the test now,
     and let the proper inverse-variance path (PR 2) retire whole-stamp `sigma_mad` wholesale — which
     resolves R1 without reviving the old machinery. If R1 proves to matter for the interim, surface
     it; **do not block on it, and do not post to Martin** — these science calls go into this fiber's
     report for Cail to fold into the eventual #741 reply.

2. **SExtractor `BACKGROUND_RMS` baseline — separate PR, closes #604.** The in-house path (THELI
   weights are an external product and slot in later as another `ME_IMAGE_PATTERN` source). Follow the
   report's change plan: `CHECKIMAGE = BACKGROUND, BACKGROUND_RMS` (config), cut the RMS stamp via
   vignetmaker's existing list-driven ME loop (config), wire one **optional** RMS input through
   `ngmix_runner`/`Vignet`, and rewrite `prepare_ngmix_weights` to build per-pixel `1/RMS²` gated on
   the mask, with `sigma_mad`/`get_noise` retained as the **fallback** when no RMS map is present.
   Validate with a before/after on the example tile (shear/s2n/flags). Effort ~4–8h.

**Quality bar:** every claim verified by running, not inferring (dev image + ngmix 2.4.0 are ready —
see Context). The red→green test must actually go red on current code and green after. **Merge and
PR-creation stay Cail's gesture** — prepare branches, commits, and a PR-ready description; surface the
science decisions in this fiber's `report.html`/outcome rather than pushing to Martin.

## Context

- **Detailed spec:** `.felt/review-ngmix-v2-pr740/weights-report.md` (change plan, file list, risks,
  the 6 science questions for Martin). This constitution is the pointer; the report is the substance.
- **R2 checkpoint:** `/tmp/pr740-wt` commit `953a52a3` restores
  `weight_map[weight_map != 0] = 1` in `prepare_ngmix_weights` and adds
  `test_weight_map_recovers_injected_inverse_variance`. Observed red before the fix:
  recovered `8.80916e11` vs truth `1e6`; observed green after the fix: recovered
  `8.80916e5` vs truth `1e6` with `rtol=0.15`. The tolerance is intentionally loose
  enough to leave R1 visible: whole-stamp `sigma_mad(gal)` still overestimates the
  compact simulated stamp's noise by about 6%, producing a ~12% inverse-variance
  shortfall. This commit fixes R2 only; it does not claim R1 is solved.
- **Current test caveat:** the targeted unit test passes in the dev image. The full
  `src/shapepipe/tests/test_ngmix.py` file currently still fails in the existing
  metacal smoke test because the container's `ngmix.fitting` has no `Fitter`
  attribute; that failure is independent of the new R2 unit path.
- **Repo:** `/automnt/n17data/cdaley/unions/shapepipe`, branch `ngmix_v2.0` (= #741 head). Base the
  work on `ngmix_v2.0` **after** the review-cleanup commit lands (it removes dead code and is currently
  staged at `/tmp/pr740-wt` commit `bd60dc8e`, pending Cail's push). Coordinate with
  [[ngmix-size-columns]] — it also edits `ngmix.py`, but a different function (the size-column writer
  vs `prepare_ngmix_weights`); low conflict, but each work-stream gets its own branch.
- **Test environment:** dev-image sandbox at `/n17data/cdaley/containers/shapepipe-dev` with **ngmix
  2.4.0** installed (the published `:develop` image ships the old `aguinot` fork v1.3.6 — wrong for
  this branch). Run via `CONTAINER_PATH=/n17data/cdaley/containers/shapepipe-dev PYTHONPATH=<worktree>/src app python3.12 -m pytest …`.
  The `--writable` sandbox needs bind mountpoints to exist (`mkdir` them); see this session's history.
- **Issue #604:** megapipe weights are 0/1 masks (Gwyn confirmed); aguinot's BACKGROUND_RMS proposal is
  the basis for PR 2. Cail confirmed the SExtractor-RMS baseline over THELI.
- **Agent:** Codex (Cail's call) — gritty, well-defined cross-file implementation.
