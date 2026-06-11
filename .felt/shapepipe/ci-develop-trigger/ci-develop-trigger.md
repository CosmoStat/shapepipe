---
id: 01KTCHWZX8GRCMJHGCRARBGFS8
name: CI silently broken on develop; install_shapepipe conda.sh lookup
status: closed
tags:
    - shapepipe
    - ci
    - bug
created-at: 2026-05-29T12:13:57.703143397+02:00
closed-at: 2026-05-29T14:09:22.682772386+02:00
outcome: Resolved by removing conda entirely (#733, merged). Enabling CI on develop exposed a multi-layer install failure in a conda env nothing deploys; rather than repair it (conda.sh lookup → --env-dev flag → unsolvable prod env → …), the conda machinery was deleted and CI now runs pytest inside the dev image. Develop is green incl. publish. Full desired state + remaining follow-ups in [[shapepipe/ci-green-on-develop]].
---

## What happened

`ci-release.yml` (the "Full Test Suite": `pytest` + the example pipeline)
only triggered on `pull_request` to `main`/`master`. Once `develop` became
the integration branch, every feature and dependabot PR targeted develop —
where CI never ran. The suite was effectively dark.

Two trigger changes landed on develop (merged directly, fast-forward):
- PR-to-develop trigger, **plus a `push: [develop]` trigger** so that
  post-merge collisions on the actual develop HEAD get tested (a
  `pull_request` run only re-fires when the PR changes, not when the base
  advances — so two independently-green PRs can break develop once both
  merge).

The very first develop run went red — including on a commit that *only*
touched the workflow trigger. That ruled out our change as the cause: the
breakage was pre-existing and merely invisible.

## The failure chain (peeled back one layer at a time)

The suite died at `./install_shapepipe`, before any test ran, on both
Linux and macOS. Each fix revealed the next failure:

**Layer 1 — conda.sh not found.** Failed instantly:
```
ERROR: Could not find /etc/profile.d/conda.sh in /opt/conda or $CONDA_PREFIX.
```
`install_shapepipe` only checked `/opt/conda` (docker convention) and
`$CONDA_PREFIX`. On runners conda lives elsewhere and `$CONDA_PREFIX` is
empty in that shell. It had passed `conda -V`, so `conda` is on PATH →
`conda info --base` resolves the real prefix. Fix keeps the old paths as
fallbacks.

**Layer 2 — wrong conda environment built.** With conda.sh found, install
ran for minutes then failed on an unsolvable env. The CI log named the env
`shapepipe` (production) and the conflict was `astropy==5.2` — pinned in
`environment.yml`, not `environment-dev.yml`. Root cause: `--develop` and
`--env-dev` are *different* flags. `--develop` sets `DEVELOP=TRUE` (adds
`.[dev]` pip extras) but leaves `ENV_DEV` false, so the build falls to the
production `environment.yml`. `--env-dev` sets `ENV_DEV=TRUE` → builds
`environment-dev.yml` as `shapepipe-dev` AND triggers `.[dev]`. The test
step's `conda activate shapepipe-dev` only ever worked with `--env-dev`.
Fix: both install steps now use `--env-dev`.

`environment-dev.yml` solves cleanly (verified locally via
`conda env create --dry-run`; floats to py3.14).

**Layer 3 — production env is genuinely broken (out of scope).**
`environment.yml` pins `astropy==5.2` (py38–py311 builds only) against
`python=3.12` → unsatisfiable. This breaks real production installs, not
just CI. Deserves its own issue; #732 doesn't touch it.

## Knock-on

**#729** (actions group, bumps `setup-miniconda`
v3→v4) hit the layer-1 failure too — confirming the action bump alone
doesn't fix the path. #729 must rebase on top of #732 once it merges before
it can go green. The smoke-test work in [[shapepipe/smoke-test-read-only]]
(#731, merged) and the lockfile group (#730, merged) landed cleanly.
