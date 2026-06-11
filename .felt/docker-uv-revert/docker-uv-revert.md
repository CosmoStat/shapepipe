---
name: 'Docker: revert skaha→python base, switch to uv lockfile'
status: active
tags:
    - shapepipe
    - docker
    - infra
created-at: 2026-04-27T11:26:45.677512058+02:00
outcome: 'PR #719 (chore: switch Dockerfile to slim Python + uv lockfile) opened and CI-green on first try (3m31s); ready for review. Drops conda double-install, makes pyproject SSOT + uv.lock the pinned manifest, switches WeightWatcher from sed-patched source build to Debian''s pre-patched 1.12+dfsg-3 package, adds binary smoke tests to deploy-image.yml.'
decisions:
    base:
        label: Base image
        rationale: Conda double-install was the actual problem; cleanest resolution is to drop conda entirely. The canfar deployment concern is satisfied as long as the slim image works on canfar.
        default: python-slim
        options:
            python-slim:
                label: python:3.12-slim with deps via uv
            skaha-astroml:
                label: 'skaha-astroml: brings conda Python; pip+conda double-install everything; image bloats 5-10 GB'
                excluded: true
                excluded_reason: Conda was the source of the bloat and double-install
            skaha-no-conda:
                label: skaha base + uv (avoid conda's python)
                excluded: true
                excluded_reason: Compromise with unclear benefit; revisit if canfar deployment specifically needs skaha
    ci:
        label: CI binary smoke test
        rationale: The original failure mode was a binary that built fine but couldn't run on canfar. Adding sextractor/weightwatcher invocation to deploy-image.yml will catch that class of regression at PR time, not deployment time.
        default: binary-smoke
        options:
            binary-smoke:
                label: Build + invoke sextractor and weightwatcher inside built image
            build-only:
                label: Build-only (current); misses runtime regressions like the weightwatcher/sextractor failure on canfar
                excluded: true
    deps:
        label: Dependency management
        rationale: We want control over when upstream changes propagate; lockfile gives that. pyproject becomes the SSOT for what shapepipe needs; uv.lock is the pinned manifest that's fully reproducible.
        default: uv-lock
        options:
            hand-pinned:
                label: Hand-pinned in Dockerfile (status quo); pyproject and Dockerfile both list deps and drift
                excluded: true
                excluded_reason: Source of repeated bugs; no upstream-pinning discipline
            poetry:
                label: Poetry + poetry.lock
                excluded: true
                excluded_reason: Same shape as uv but ~10x slower; uv is the modern default
            uv-lock:
                label: uv + pyproject + uv.lock; uv sync --frozen in Dockerfile
    modernize:
        label: Modernize package versions
        rationale: 'We determined which versions MUST stay pinned: only ngmix (pinned to a stable_version fork branch — replacement is tracked separately). Everything else can move to current latest because uv resolved cleanly and CI smoke test still passes (3m42s). If a real pipeline run on canfar surfaces a numpy-2 / pandas-3 break, the fix is a targeted constraint + uv lock, not a wholesale revert.'
        default: stay-current
        options:
            stay-conservative:
                label: Keep pre-v2 minimums (numpy 1.26, astropy 6.1, pandas 2.2); only bump when forced
                excluded: true
                excluded_reason: Drift between pyproject signal and lockfile reality; loses the chance to surface numpy-2/pandas-3 incompatibilities at PR time when CI is fast
            stay-current:
                label: Bump pyproject minimums to current major versions (numpy 2, astropy 7, pandas 3, galsim 2.8, mpi4py 4.1, etc.); pin ngmix to its stable_version fork branch
insights:
    ci-fast:
        claim: 'First CI run on PR #719 went green in 3m31s. uv installed 238 packages in 322ms — everything resolved to prebuilt wheels, no source compilation of galsim/mpi4py/python-pysap/etc. Massive speedup vs. previous build.'
    failure-mode:
        claim: The triggering failure was a binary (weightwatcher / sextractor) that built fine but couldn't run on canfar. Build-only CI doesn't catch that class; binary-version CI does.
    patches-explained:
        claim: WeightWatcher 1.12 (2014) source uses pre-GCC-10 'common symbol' globals (prefstruct prefs; etc. in headers). GCC 10 (Apr 2020) flipped default to -fno-common, breaking the build. Old Dockerfile fixed it inline with sed converting headers to extern + adding single .c definition. Debian's weightwatcher 1.12+dfsg-3 package applies an equivalent patch upstream of the build, which is why apt install just works.
    skaha-uv-incompatible:
        claim: skaha-astroml's value comes from conda's pre-built scientific stack; layering uv on top defeats that purpose, so "skaha + uv (no conda)" is incoherent rather than just hard.
    weightwatcher-apt:
        claim: Debian bookworm packages WeightWatcher 1.12 directly — no need to build from source with the patched headers the v2 Dockerfile uses.
---

## Where the work landed

[CosmoStat/shapepipe#719](https://github.com/CosmoStat/shapepipe/pull/719) — draft PR.

Commit: `9af10257 chore: switch Dockerfile to slim Python + uv lockfile`. Four files:

- `Dockerfile` — `python:3.12-slim-bookworm` base, apt for astromatic binaries (psfex, source-extractor, weightwatcher), uv for Python deps.
- `pyproject.toml` — bumped minimums to current versions; new `jupyter` extra holds ipython/jupyterlab/snakemake (previously inlined in Dockerfile).
- `uv.lock` — 286 packages resolved by `uv lock`.
- `.github/workflows/deploy-image.yml` — split `Test` into `Test — binaries` and `Test — shapepipe entry point`.

## Lockfile workflow notes

```bash
uv lock                # regenerate after editing dependencies
uv add 'pkg>=1.2'      # adds to pyproject AND uv.lock
uv lock --upgrade      # bump pinned versions deliberately
uv sync --frozen       # in Dockerfile; fails if pyproject/uv.lock drift
```

The `--frozen` flag is the discipline mechanism: a stale lockfile cannot ship.

## Followups

- Watch CI on #719. The slim-base apt list is conjectural — galsim/mpi4py/python-pysap pull a lot of system deps and we may need to add more (`libatlas-base-dev`, `libblas-dev`, etc).
- If CI needs anything beyond what's in the apt block, that's worth noting for next time.
- After this lands, PRs #708 and #714 may need a small rebase.
- Optional: separate `Dockerfile.canfar` building on skaha if there's a concrete deployment reason. Currently conjectural — floated as a possibility, but slim should work on canfar.

## Connections

- [[shapepipe]] — root
