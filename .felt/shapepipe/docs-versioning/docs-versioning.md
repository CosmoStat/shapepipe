---
name: Versioned docs site with a version switcher
tags:
    - shapepipe
    - docs
    - ci
created-at: 2026-05-31T20:39:06.459699229+02:00
outcome: 'PR #738 (off develop): docs site now versioned — stable@root (non-destructive to existing URLs), develop@/latest/, tags@/vX.Y.Z/, switcher.json drives the dropdown. Fixes that the site published only from master (~78 commits behind) so container.md/canfar.md were never published. CI green incl. real dev-image build + PR-preview artifact.'
---

## Why

`cd.yml` published the Sphinx site only on push to `master`, which sits ~78
commits behind `develop`. So the public docs described a stale release while
everyone runs `develop`, and six pages written on `develop` — `container.md`,
`canfar.md`, four others — had **never been published** (`container.html` →
404). Surfaced when checking where the README's container link would land.

## Design

Versioned, not just "publish from develop." gh-pages layout is **additive**:
`master`→ site **root** (stable; keeps every existing `/installation.html`-style
URL working), `develop`→ `/latest/`, tags `v*`→ `/<tag>/`. A root `switcher.json`
(tracked in `docs/`) drives the pydata/`sphinx-book-theme` dropdown;
`DOCS_VERSION` (set by the workflow per ref) tells `conf.py` which entry to
highlight. Each deploy uses `keep_files: true`, so it only writes its own slice
and never clobbers siblings — every step is non-destructive, including the first
`develop` run (adds `/latest/` + `switcher.json`, leaves root untouched).

**Per-ref builds, not `sphinx-multiversion`.** autodoc imports `shapepipe`, so an
old tag must build against *its own* deps. Multiversion uses one shared env
across all refs → breaks autodoc on old tags. Building each ref in its own dev
image at its own time sidesteps that entirely. PRs touching docs build + upload
the HTML as an artifact but don't deploy → broken builds caught pre-merge,
reviewers get a preview. `conf.py` version now comes from package metadata
(killed the hand-maintained `1.0.1`, actual is `1.1.0`).

## The recurring bit-rot pattern

Same shape as [[shapepipe/mpi-hybrid]]: a path the repo's CI never exercises
silently rots, and *enabling* it surfaces the breakage. Here, the docs build
only ever ran on `master` (old flat layout), so `sphinx-apidoc … shapepipe`
went stale when the package moved to src-layout (`src/shapepipe`); the first
`develop` build failed `/app/shapepipe is not a directory`. Fixed to
`src/shapepipe`. The lesson keeps repeating — unexercised paths (MPI, docs
deploy) drift undetected; the fix is to *run them in CI*, which this PR now does
for docs on every push and every docs-touching PR.

## Follow-ups (not in PR)

Seed `switcher.json` with the first real release entry when the next version is
tagged; one-time `master` release to refresh root with current docs. Policy
nod: the site now publishes `develop` as `latest`, so in-development docs are
publicly visible — right for an active pipeline, but a conscious shift.
