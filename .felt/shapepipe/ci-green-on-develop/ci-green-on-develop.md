---
id: 01KTCHWZX9Q1MG2FB20N5Y52TD
name: 'ShapePipe CI: green & trustworthy test suite on develop'
status: open
tags:
    - shapepipe
    - ci
    - constitution
created-at: 2026-05-29T12:34:31.800189806+02:00
outcome: 'Near-realized. Ground truth 2026-06-11: develop green incl. image publish; #708 closed (test scaffolding landed on develop: 5 structural test files + property/synthetic-FITS tiers); candide PBS scripts modernized to SLURM+apptainer (#737); CLAUDE.md already a solid onboarding doc (minor gaps: CONTRIBUTING pointer, first-PR walkthrough); docs deploy infra correct post-#738 — /latest/ live, switcher published, ROOT still serves 2022 v1.0.1 until the next master push (the one remaining trigger). Conda survives only in out-of-scope CANFAR/CC-IN2P3 scripts (job_sp_canfar.bash + init_run_exclusive_canfar.sh have ACTIVE conda paths; cc_{mpi,smp}.sh use ccenv anaconda) — needs its own cluster-aware pass to fully realize ''no conda anywhere''.'
---

## Desired State

**There is no conda anywhere in the repo.** The Dockerfile is the one
dependency story: `python:3.12-slim-bookworm`, astromatic binaries
(`psfex`, `source-extractor`, `weightwatcher`) + system libs via `apt`,
Python deps via `uv` from `pyproject.toml` + `uv.lock` (the SSOT), shapepipe
via `uv pip install -e .`. Python is pinned to 3.12 by the base image and
`requires-python>=3.12`.

**CI tests the artifact it ships.** The dev image (runtime + all extras,
including `test`) is built, then `pytest` runs *inside it* — alongside the
example-pipeline smoke and binary checks already in `deploy-image.yml`. A
green check means the actual container works, not a parallel conda env that
nothing deploys.

**The workflow separates concerns:**
- **Test** — build dev image, run `pytest` + smoke (example pipeline,
  binaries) in it. Runs on PRs to develop/main/master *and* pushes to
  develop. No registry push on PRs.
- **Publish** — push runtime + dev images to ghcr. Only on push to
  develop/main, and only after Test passes.

**Docs build from the image too** — `cd.yml` (deploy API docs to gh-pages on
master) and whatever survives of `doc-tests.yml` use the dev image's `doc`
extra (sphinx, myst, …) instead of conda + `install_shapepipe`. (Watch for
`pandoc`, which was `conda install`ed; add it to the dev image via apt if the
build needs it.)

**Done conditions**
- `install_shapepipe`, `environment.yml`, `environment-dev.yml` deleted.
- No `setup-miniconda` / `conda` / `install_shapepipe` references in build
  machinery, CI, or install docs. (Cluster job scripts under `example/pbs/`
  and `scripts/sh/*.bash` still hardcode personal conda env paths — those are
  operational and cluster-specific; fenced out below, not a blocker.)
- `pytest` green inside the dev image on develop (Linux).
- Docs still build and deploy.
- `CONTRIBUTING.md` and `docs/source/installation.md` document the
  container/uv workflow, not conda.

**Scope fence**
- **Linux only.** Deployment is Linux apptainer/SIF; the dev container is
  also the Mac developer story. The macOS CI leg is dropped, not ported to
  Homebrew.
- Not expanding the test suite (that's #708 scaffolding) and not fixing
  scientific-module bugs `pytest` may reveal — if one blocks green, xfail/skip
  with a tracked issue rather than a deep fix here.
- **Cluster job scripts** (`example/pbs/*.sh`, `scripts/sh/job_sp_canfar*.bash`,
  `init_run_exclusive_canfar.sh`, …) hardcode `~/.conda/envs/shapepipe` paths.
  Porting them to the apptainer/container workflow is a separate, cluster-aware
  pass — touching Martin's operational scripts blind would break live runs.
- The production `environment.yml` `astropy==5.2` bug is **moot** — the file
  is being deleted, not fixed.

## Context

How we got here, and the diagnostic layers that justified deletion over
repair, live in [[shapepipe/ci-develop-trigger]]. This fiber is the desired
state; that one is the autopsy.

The conda-free path is already proven by the Dockerfile — read it first; it
*is* the spec for what CI and docs should do. Key facts:
- Astromatic tools are Debian packages (`apt-get install psfex
  source-extractor weightwatcher`); that was conda's whole reason for being.
- `pyproject.toml` has a `dev` extra = `doc,jupyter,lint,release,test,fitsio`;
  the dev image pre-installs it. `[project.scripts]` defines the console
  entry points.
- `deploy-image.yml` already builds both targets, runs the example pipeline
  read-only (#731) and binary smokes; it only lacks the real `pytest` run
  (currently just `pytest --version`) and a test/publish split. It triggers
  on every branch `push` and publishes unconditionally — the split fixes
  that.

Files still carrying conda to remove/convert: `ci-release.yml` (delete),
`cd.yml` + `doc-tests.yml` (convert to image-based docs build),
`CONTRIBUTING.md`, `docs/source/installation.md`.

In-flight PRs affected: **#732** (conda.sh + `--env-dev` fix to
`ci-release.yml`) is now superseded — close it, we're deleting that file.
**#729** (dependabot actions group) bumps actions in workflows we're
editing/deleting — rebase it after this lands so it only touches survivors.

## Skills

`felt`. Verify empirically: the real GH workflow run and `docker build` +
`docker run … pytest` locally are ground truth, not inference.

## Evidence

```bash
# No conda left anywhere?
grep -rIn -e conda -e install_shapepipe -e setup-miniconda -e environment-dev \
  --exclude-dir=.git --exclude-dir=.felt .

# Does the dev image build and pass the suite locally?
docker build --target dev -t sp-dev . && \
  docker run --rm sp-dev pytest

# CI signal on develop
gh run list --branch develop --limit 6 \
  --json workflowName,event,status,conclusion \
  --jq '.[] | "\(.workflowName) [\(.event)] \(.status)/\(.conclusion)"'
```

Done = the grep is clean, the suite is green inside the image on develop,
and docs still deploy.

## Open Questions

- **`doc-tests.yml` — convert or delete?** It's `workflow_dispatch`-only and
  largely duplicates `cd.yml`'s build. Likely fold into a single image-based
  docs job rather than maintain two.
- **Does the docs build still need `pandoc`?** It was `conda install`ed in
  `cd.yml`. If sphinx still needs it, add `pandoc` to the dev image (apt);
  if it was vestigial, drop it.
- **Bare-metal cluster installs** — anyone who ran `install_shapepipe`
  directly (no container) loses it. Is the container/apptainer path
  sufficient for all cluster users, or does a documented non-container
  install need to survive in `installation.md`?
