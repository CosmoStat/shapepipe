---
id: 01KTCHWZXE6Z6MG1P21F804RZ2
name: 'ShapePipe cleanup: remove obsolete rho-stats/stile; modernize candide job scripts'
status: closed
tags:
    - shapepipe
    - cleanup
    - constitution
created-at: 2026-05-30T21:45:50.977369486+02:00
closed-at: 2026-05-30T22:30:55.745891097+02:00
outcome: |-
    Mixed. PR #737 (candide scripts) good and open; PR #736 (rho-stats) CLOSED as
    an over-cut. The constitution's premise was stale: it treated mccd_plots_runner
    as "the rho-stats path," but the authorized rho-stats removal (Martin #657
    "option 1, delete") had ALREADY shipped in #715 (merged 2026-04-23). By develop,
    mccd_plots_runner/mccd_plot_utilities held NO rho/stile/treecorr — pure
    meanshapes (plot_meanshapes: focal-plane ellipticities/sizes/residuals + 1D
    histograms). On #715 Martin EXPLICITLY said keep meanshapes ("very useful...
    maybe Fabian's 2D plots in the catalogue paper"). #736 deleted exactly that path
    → net-zero rho cleanup, entirely over-reach against Martin's instruction. Closed
    with explanation (2026-05-31), not trimmed — nothing left to salvage.
    (2) PR #737 stands: candide_smp.sh / candide_mpi.sh run via apptainer + runtime
    image, no conda; SMP verified end-to-end on c03 (0 errors), MPI hybrid pattern
    written but needs a real allocation to verify. canfar + ccin2p3 left untouched.
    Scope findings that held up: `stile` already vestigial; `random_cat` correctly
    KEPT (general LSS random generator, not rho-stats). LESSON: check what prior
    merged PRs already did before acting on a cleanup constitution — verify the
    premise against current `develop`, not the constitution's snapshot.
shuttle:
    enabled: true
    kind: oneshot
    host: c03
    project_dir: /automnt/n17data/cdaley/unions/shapepipe
    agent: claude-opus
    session:
        id: 7c02b521-17bd-44dc-ad10-709520d6d2bf
        agent: claude-opus
        dispatched_at: 2026-05-30T20:29:24.913621662Z
tempered: true
---

## Desired State

Two separate, reviewable PRs against `develop` — **opened, not merged.** This is
cleanup someone else (Martin) should sign off on; the worker's job is to land
the change on a branch and open the PR, not to merge it.

### Deliverable 1 — remove the obsolete rho-stats / stile path (PR for Martin)

The in-ShapePipe rho-statistics path is dead code. Martin confirmed on #657
(2026-04-22): *"the rho stats within shapepipe are obsolete, and have been
superseded by sp_validation / cosmo_val, which uses shear_psf_leakage."* The
stack division is: ShapePipe **produces** catalogues, `sp_validation`
**validates**. So this path should be deleted, not migrated or preserved
"just in case."

Done state:
- The rho-stats modules are gone: `modules/mccd_plots_runner.py`,
  `modules/mccd_package/mccd_plot_utilities.py`, and the random-catalogue
  null-test path (`modules/random_cat_package/`, `modules/random_cat_runner.py`)
  — **after confirming** `random_cat` is only used for rho-stats nulls (it
  appears in `example/cfis/config_Rc.ini`; check it isn't load-bearing for any
  non-rho clustering use before deleting).
- Whatever pulls `stile` is removed (`stile` is *not* a declared dependency in
  `pyproject.toml`/`uv.lock` and isn't directly imported in those modules — so
  first find what actually references it; it may already be vestigial. Remove
  any real reference; don't invent one).
- The **9 example configs** that name `mccd_plots_runner` in their `MODULE =`
  line — including `example/cfis_simu/config_tile_Sx_exp_mccd.ini` (the sim
  config) — no longer reference it: drop it from the `MODULE` lists and remove
  the now-orphaned `[MCCD_PLOTS_RUNNER]` sections. Same for `random_cat_runner`
  / `[RANDOM_CAT_RUNNER]` if random_cat goes. Per Martin's note, assume the
  `RHO_STATS_STYLE = UNIONS` path is also not load-bearing.
- `docs/source/random_cat.md` (and any rho-stats doc) is removed/updated to match.
- The in-image test suite stays green — the structural tier (config-parse,
  imports, runner-metadata) will flag a broken config or dangling import.

### Deliverable 2 — modernize candide job scripts (PR, tested here)

`example/pbs/candide_smp.sh` and `candide_mpi.sh` activate a personal conda env
(`module load intelpython/3; source activate $HOME/.conda/envs/shapepipe`) and
call `$SPENV/bin/shapepipe_run`. Convert them to run the pipeline through the
**container** (apptainer, the published `ghcr.io/cosmostat/shapepipe` image)
instead — no conda. **This host is c03 = candide, so test the converted script
here** (at minimum, the apptainer invocation runs `shapepipe_run` on the example
config; a full PBS submission if practical).

The **canfar** scripts (`scripts/sh/job_sp_canfar*.bash`,
`init_run_exclusive_canfar.sh`, etc.) are **out of scope** — they target a
different cluster we can't verify here. Leave them untouched and say so in the PR
description.

## Context

- **Why rho-stats is dead:** see the memory note (Martin #657 comment) and the
  stack division — `sp_validation`/`cosmo_val` owns PSF-systematics validation
  now via `shear_psf_leakage`. Don't preserve the UNIONS rho path out of caution.
- **Files in the rho path:** `modules/mccd_plots_runner.py`,
  `modules/mccd_package/mccd_plot_utilities.py`, `modules/random_cat_package/`,
  `modules/random_cat_runner.py`.
- **Configs to edit** (grep `mccd_plots_runner` / `random_cat_runner` in
  `example/`): `config_Rc.ini`, `config_MsPl_psfex.ini`,
  `config_valjoint_Pl_mccd.ini`, `config_exp_mccd.ini`, `config_MsPl_stars.ini`,
  `config_Pl_mccd.ini`, `config_MsPl_mccd.ini`, `cfis/defunct/…`, and
  `example/cfis_simu/config_tile_Sx_exp_mccd.ini`.
- **CI:** the only gate is `pytest` inside the dev image (`deploy-image.yml`);
  run it in the image to verify (`docker run --rm ghcr.io/cosmostat/shapepipe:develop pytest`).
- **Container for the candide scripts:** the runtime image
  (`:develop-runtime`) is the slim target for batch jobs; see
  `docs/source/container.md` and `docs/source/installation.md` for the apptainer
  invocation pattern.

## Skills

`felt`. Verify empirically — run the in-image suite for D1; run the converted
candide script via apptainer on c03 for D2.

## Evidence

```bash
# D1: nothing references the removed rho path
grep -rIn "mccd_plots_runner\|random_cat_runner\|MCCD_PLOTS_RUNNER" src/ example/
# D1: suite still green in the image
docker run --rm ghcr.io/cosmostat/shapepipe:develop pytest -rX
# D2: the converted candide script runs the pipeline via the container, no conda
grep -rn "conda\|source activate\|SPENV" example/pbs/candide_*.sh   # → none
```

Done = two PRs open against develop (rho-stats removal; candide scripts), CI
green on each, neither merged, canfar untouched-and-noted.

## Scope

- **In:** the rho-stats/stile code + its config references + docs; the candide
  PBS scripts.
- **Out:** merging (open PRs only — Martin reviews); the canfar scripts (note,
  don't touch); the broader test-suite work (separate, done); any
  scientific-algorithm change.

## Resolution

**`random_cat` is NOT rho-stats-only — kept.** Empirically: `random_cat_package`
generates a random catalogue of points within the survey mask (healpix output,
`config_Rc.ini`, authored by Martin). That is a general clustering / LSS tool
(Landy-Szalay-style randoms), independent of PSF systematics. Rho statistics are
auto/cross-correlations of PSF-ellipticity residuals at *star* positions and need
no random catalogue. So the conditional in Deliverable 1 ("remove after confirming
random_cat is rho-only") failed its precondition → random_cat stays. Removing it
would overreach what Martin called obsolete; flagged in PR #736 for a follow-up if
he does want it gone.

**`stile` was already vestigial.** Zero references in `src/`, `example/`, `docs/`,
`pyproject.toml`, `uv.lock` — never a declared dep, never imported. Nothing to
remove; the absence is noted in PR #736 so it isn't read as an oversight.

**The "rho-stats path" was really PSF-diagnostic plotting.** `mccd_plot_utilities`
imports only matplotlib/mccd/numpy/astropy — it computes mean shapes and ellipticity
histograms, no `treecorr`/`stile` correlation functions. The docstring's "rho
statistics plot" was aspirational. Maps cleanly to Martin's "rho stats within
shapepipe" = the PSF-leakage diagnostics now owned by `shear_psf_leakage`.

**MPI verification gap.** `candide_smp.sh` is verified end-to-end through the
container on c03. `candide_mpi.sh` uses the documented hybrid Apptainer pattern
(host `mpiexec` launches container ranks) but can't be verified on a login node —
`mpiexec` inside the container hangs without a PMIx allocation. ABI note for the
reviewer: image ships OpenMPI 4.1.4, candide modules are now OpenMPI 5.0.x.
