---
name: ShapePipe — project knowledge & active threads
tags:
    - shapepipe
created-at: 2026-04-27T11:26:38.71538657+02:00
outcome: 'Root of ShapePipe''s felt store: the stack division, repo conventions, and the why behind in-flight infra/cleanup threads.'
---

This is the root of ShapePipe's felt store — shared notes on architecture
decisions, conventions, and in-flight work, for the team and AI agents alike.
ShapePipe is the UNIONS galaxy shape-measurement pipeline; `CLAUDE.md` covers the
build / container / CI overview, and the fibers here carry the *why*. Start here,
then follow the links.

## Stack division

ShapePipe **produces** shear catalogues; `sp_validation` / `cosmo_val`
**consume** and validate them; `cs_util` holds code shared across both. A concern
about *validating* catalogues belongs downstream, not in ShapePipe.

## Conventions specific to this repo

- **Rho-statistics are obsolete inside ShapePipe.** PSF-systematics validation
  moved downstream to `sp_validation` / `cosmo_val` (via `shear_psf_leakage`);
  the stile/treecorr rho code was removed in #715. But the **meanshapes /
  ellipticity focal-plane plots** (`mccd_plots_runner`) are *deliberately kept* —
  they are a general PSF/star-catalogue diagnostic, not rho-stats, and feed
  catalogue-paper figures. Don't delete that path along with rho-stats; see
  [[shapepipe/cleanup-rhostats-jobscripts]] for where the boundary actually sits.
- Run the pipeline through the container; use `python3.12` explicitly inside it.
- **ngmix** is pinned to a fork branch until fixes land upstream — don't bump
  that dependency line. [[ngmix-update]] tracks the path back to upstream.

## Active threads

- **[[shapepipe/ci-green-on-develop]]** / **[[shapepipe/test-suite]]** — a
  tiered, in-image test suite and trustworthy CI on `develop`.
- **[[docker-uv-revert]]** — slim Python base + uv lockfile, dropping conda.
- **[[shapepipe/mpi-hybrid]]** — running hybrid MPI through the container on candide.
- **[[ngmix-update]]** — replacing the pinned ngmix fork with upstream.
