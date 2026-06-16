---
id: 01KTCHX00NMQ1VDPGXRYJ6RGZR
name: ShapePipe maintenance & PRs
tags:
    - shapepipe
    - portolan
created-at: 2026-04-27T11:26:38.71538657+02:00
outcome: 'Root: collaboration with Martin on ShapePipe — PRs, infra, future ngmix and Fabian work'
---

ShapePipe is the UNIONS shape-measurement pipeline. I'm not the primary
maintainer (that's Martin Kilbinger); my role is collaborator helping
clean up infra, surface bugs, and keep the merge queue moving while
Martin focuses on science threads.

## Working agreement with Martin

Surfaced over a 2026-04-27 walking conversation. Captured in
[[shapepipe/prs-in-flight]] and the per-thread fibers below.

- I review and patch his PRs; he reviews mine. Bugs found during review
  go to a dedicated PR rather than getting bundled into his feature
  branch (per `feedback_separate_infra_prs`).
- v2.0 was merged fast (it was ready). The skaha base it brought in is
  the active source of pain → see [[shapepipe/docker-uv-revert]].
- I file the issues; Claude usually drafts the PRs in my voice.
  Disclosure on Claude-only review per
  `feedback_claude_only_review_disclosure`.

## Active threads (refreshed 2026-06-10)

- **PR #737** (MPI on candide + containerized SLURM scripts) — **MERGED 2026-06-10** (rebased onto develop, both CI runs green).
- **PR #738** (versioned docs + switcher) — **MERGED 2026-06-10** after independent review; sfarrens offered post-hoc comments. Stable root refreshes on next master push. See [[docs-versioning]].
- **PR #739** (machine-specific cluster docs tree) — awaiting Martin's review; check rebase state after the #737/#738 merges. See [[docs-cluster-tree]].
- **PR #741** (ngmix v2.0, CI mirror of #740) — Martin left 10 inline comments Jun 5, no verdict yet; shear recovery verified unbiased to ~1e-4 in m. See [[ngmix-weights-ivar]] (PR-ready regression fix + #604 ivar plan) and [[ngmix-size-columns]] (honest r50 spec).
- **Martin's PRs #704 (contributors) & #699 (coverage mask)** — Cail's review requested; both conflicting, need Martin's rebase first.
- **[[ci-green-on-develop]]** — conda fully removed (#733 merged); remaining follow-ups tracked there.

## Earlier threads (superseded)

- [[shapepipe/docker-uv-revert]] / [[shapepipe/prs-in-flight]] — the #719-era conda removal; landed via #733, current state in [[ci-green-on-develop]].
- **[[shapepipe/ngmix-update]]** — became the #740/#741 ngmix v2.0 review thread.
- **[[shapepipe/fabian-coord-bug]]** — port Fabian's 1-line coord propagation fix; still pending his image-sim code reaching github.

## Conventions specific to this repo

- Container runs through `app` (apptainer wrapper); use `python3.12`
  inside the shapepipe container (see `reference_containers`).
- ShapePipe produces; `sp_validation` consumes; `cs_util` is shared (see
  `project_stack_division`).
- Rho stats are obsolete here — sp_validation/cosmo_val took over (see
  `project_rho_stats_obsolete`).
- Royal "we" in PR/issue voice; specific findings attributed to Claude
  by name (see `feedback_writing_voice_on_cails_behalf`).
