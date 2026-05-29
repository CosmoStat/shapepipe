---
name: Triage open dependabot uv.lock PRs
status: closed
tags:
    - shapepipe
    - constitution
    - pr
    - dependabot
created-at: 2026-05-27T11:57:31.076323766+02:00
closed-at: 2026-05-28T15:27:23.026Z
outcome: All 6 PRs merged in-conversation (727 idna, 722 mistune, 724 gitpython, 720 jupyterlab + pyproject floor, 721 jupyter-server, 726 urllib3). 19 of 20 Dependabot alerts cleared. Last one (sqlitedict) dismissed as tolerable_risk; see [[shapepipe/sqlitedict-pickle-smell]]. Shuttle daemon didn't pick up the dispatch — root cause unknown, worth a separate look.
shuttle:
    enabled: true
    kind: oneshot
    project_dir: .
    agent: pi-gpt-5.4
tempered: true
---

## Why this exists

`uv.lock` landed in `develop` on 2026-05-05 via PR #719. Before that there was no lockfile and Dependabot had nothing to pin against. After that landing, Dependabot's security-updates channel (which fires automatically without a `dependabot.yml` config) started opening one PR per GHSA-touched package in the lockfile. Six have accumulated; none merged. CI is green on the ones spot-checked. Maintainers have write access and merge their own PRs (e.g. #718), so these can be merged directly.

This fiber's job is the one-time backlog clearance. The longer-term policy (a `.github/dependabot.yml` that groups bumps, schedules them, scopes ecosystems) is **out of scope** — that's being decided in parent conversation.

## Desired State

All six PRs (#720, #721, #722, #724, #726, #727) closed — each either squash-merged into `develop`, or closed with a clear reason posted as a PR comment.

`gh pr list --repo CosmoStat/shapepipe --author "app/dependabot" --state open` returns empty when this is done.

Each merged PR's commit appears on `origin/develop` as a single squash commit with a sensible message.

## Scope (what to leave alone)

- **Do not write `.github/dependabot.yml`.** Policy is being decided elsewhere; writing it from here would race with that conversation.
- **Do not modify `pyproject.toml`.** If a bump would require widening or narrowing a direct-dep ceiling, leave a comment on the PR explaining the gap and move on. Don't push edits.
- **Do not touch non-dependabot PRs or issues.** There's other backlog (issues #709, #711, #712 from April; older issues from Martin/Lucie/Axel) — that's a separate sweep.
- **Do not edit the lockfile manually.** If a Dependabot PR has gone stale (lockfile conflict against develop), trigger a rebase via `@dependabot rebase` and move on; don't resolve conflicts by hand.
- **Don't merge across pyproject ceilings.** If `pyproject.toml` has a direct constraint like `package<X` and the bump would cross it, that's not a clean merge — comment and skip.

## The six PRs

| # | Package | Bump | Notes |
|---|---------|------|-------|
| 720 | jupyterlab | 4.5.6 → 4.5.7 | Direct dep (`jupyterlab>=4.3` in pyproject's `jupyter` extra), patch bump |
| 721 | jupyter-server | 2.17.0 → 2.18.0 | Transitive |
| 722 | mistune | 3.2.0 → 3.2.1 | Transitive, patch bump |
| 724 | gitpython | 3.1.47 → 3.1.50 | Transitive, replaced earlier #723 (3.1.49) |
| 726 | urllib3 | 2.6.3 → 2.7.0 | Transitive |
| 727 | idna | 3.13 → 3.15 | Transitive, CVE-2026-45409 (medium severity DoS) |

All sampled PRs are lockfile-only (3+/3- in `uv.lock`). Spot-checked CI was green on #727 and #720; assume the others are too but verify.

## Per-PR workflow

For each PR, in this order (lowest-risk-first: idna's CVE is the cleanest justification, lockfile-only patches before minors):

```
727 → 722 → 724 → 720 → 721 → 726
```

1. **Inspect.** `gh pr view <N> --repo CosmoStat/shapepipe --json files,additions,deletions,mergeable,mergeStateStatus`. Confirm `files` is only `uv.lock`, additions+deletions are small (~3-30).
2. **Check CI.** `gh pr checks <N> --repo CosmoStat/shapepipe`. Required: `build-and-push-image` passes. If failing, read logs (`gh run view <run-id> --log-failed`) to understand whether it's a real regression or flaky.
3. **Check freshness.** If `mergeStateStatus` is `BEHIND` or shows conflicts, comment `@dependabot rebase` and move to the next PR; loop back later.
4. **Verify no pyproject ceiling violation.** `grep -E "<package_name>" pyproject.toml` — if there's a direct upper-bound that the bump crosses, stop and comment.
5. **Squash-merge.** `gh pr merge <N> --repo CosmoStat/shapepipe --squash`. Default commit title is fine (`Bump <pkg> from X to Y (#N)`). The repo doesn't auto-delete merged branches (`delete_branch_on_merge: false`) — leave that as-is, Dependabot handles its own branches.
6. **Brief pause.** After each merge, the other dependabot branches go stale on `uv.lock`. Dependabot auto-rebases on `develop` updates, but it takes a moment. Move to the next PR; if `mergeStateStatus` shows behind, `@dependabot rebase` and circle back.

## Evidence (how to know you're done)

- `gh pr list --repo CosmoStat/shapepipe --author "app/dependabot" --state open` returns `no open pull requests`.
- `git log origin/develop --oneline -10` shows the squash-merged bump commits.
- Any PR closed without merge has a comment explaining why (e.g. "deferring, pyproject ceiling at `<X`").

## Open Questions

None for the worker — kick everything that needs human judgment back to the PR as a comment and skip that PR. The parent conversation will handle anything left after the loop.

## Handoff (when exiting)

Update this fiber's `outcome:` to reflect what landed (e.g. "5 merged, 1 deferred for pyproject ceiling on `<pkg>`"), set `status: closed`, leave `tempered` unset (the maintainer will accept once verified). Append a `felt history` event with the merge SHAs and any deferred-PR rationales.

## Skills

Activate `felt` (you'll touch the fiber). The `shuttle` skill is already loaded by virtue of being the dispatched worker.
