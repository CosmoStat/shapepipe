---
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

## Active threads

- **[[shapepipe/docker-uv-revert]]** — slim Python + uv lockfile, drop conda. PR #719 (draft).
- **[[shapepipe/prs-in-flight]]** — tracking #708 (testing scaffold), #714 (develop bugs), #719 (this one).

## Future work

- **[[shapepipe/ngmix-update]]** — replace Axel's stable_version fork
  with upstream ngmix; reconcile with Lucy's wrapper.
- **[[shapepipe/fabian-coord-bug]]** — port Fabian's 1-line coord
  propagation fix; first need his image-sim code on github.

## Conventions specific to this repo

- Container runs through `app` (apptainer wrapper); use `python3.12`
  inside the shapepipe container (see `reference_containers`).
- ShapePipe produces; `sp_validation` consumes; `cs_util` is shared (see
  `project_stack_division`).
- Rho stats are obsolete here — sp_validation/cosmo_val took over (see
  `project_rho_stats_obsolete`).
- Royal "we" in PR/issue voice; specific findings attributed to Claude
  by name (see `feedback_writing_voice_on_cails_behalf`).
