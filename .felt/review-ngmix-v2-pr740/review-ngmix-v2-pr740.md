---
name: 'Review + work: ngmix v2.0 (PR #740)'
status: closed
tags:
    - constitution
    - shapepipe
    - ngmix
    - review
created-at: 2026-06-01T11:50:17.967872934+02:00
closed-at: 2026-06-03T16:28:05.521094506+02:00
outcome: 'Review of PR #740 (ngmix v2.0) delivered as a coherent two-part review on CI-mirror #741 (part 1: fixed+tested items; part 2: empirical verification + 11 line-anchored findings + methodology questions for Martin/Lucy). Verified empirically on candide against real ngmix 2.4.0: API correct, metacal recovers m=+2e-4 (consistent w/ zero), centroid fix benign (its bias lives in old ngmix-1.x path, can''t repro on current code — expected). Key findings: weight-norm change, *_psfo now reconvolved-PSF + per-type collapse, zero-pixel guard any->all, CHECK_EXISTING_DIR resume dropped, r50/T mislabel, runner decorator contracts unupdated for 4 runners. Martin to go through and merge; offer to push cut-and-dry fixes stands. Merge is Cail''s/Martin''s gesture.'
horizon: now
shuttle:
    enabled: true
    kind: oneshot
    host: candide
    project_dir: /automnt/n17data/cdaley/unions/shapepipe
    agent: claude-opus
---

# Review + work: ngmix v2.0 (PR #740)

Martin requested Cail's review on **[CosmoStat/shapepipe #740 "Ngmix v2.0"](https://github.com/CosmoStat/shapepipe/pull/740)** (head `ngmix_v2.0` → base `develop`, +3894/−3410 across 67 files). This is the PR that **realizes [[ngmix-update]]**: replace Axel's `stable_version` ngmix fork with **upstream ngmix 2.4.0** and adopt **@lbaumo's (Lucy Baumont's) new ngmix classes + interface**, reconciling the cleaned-up wrapper from her visit. Track it in [[prs-in-flight]]; it sits in the broader Martin collaboration under [[shapepipe]].

## Desired State

PR #740 has a **thorough, scientifically-grounded review**, and the collaboration moves forward:

1. **A structured review exists**, worked against Martin's reviewer checklist (targets `develop` ✓, labels, clear description, "closes #" links, code/doc style, docs updated, **CI passing**, all changed files read with comments) — *and* against a **science-quality bar** that the checklist doesn't capture (below).
2. **Comments are posted to GitHub** in Cail's voice (signed "— Claude on behalf of Cail" per the authorship convention, or Cail's own voice as he chooses on the call) — line-level where it's a specific code point, PR-level for the broader read.
3. Where Cail decides, **fixes/suggestions are pushed** (as PR review suggestions or a branch off `ngmix_v2.0`) or **requested** from Martin/Lucy. **Merge/approve is Cail's gesture, never the worker's.**

### The science-quality bar (what the checklist misses)

The diff is dominated by `src/shapepipe/modules/ngmix_package/ngmix.py` (~1252 lines changed) — the module overhaul. The review must actually reason about the science, not just the diff mechanics:

- **ngmix 2.4.0 API correctness.** The 1.x→2.x ngmix API changed substantially (observation/fitter/bootstrapper interfaces). Does the new code use the 2.4.0 API correctly — priors, PSF fitting, the metacal/bootstrapping path, guesses? @lbaumo's new classes: are the abstractions sound, and does the interface match how shapepipe calls them (`ngmix_runner.py`)?
- **The centroid-bias fix.** A large part of this PR is centroid-bias work (`scripts/validation/centroid/centroid_bias{,_v2}.py`, the bug/fix/v2 configs, `test_centroid_shift.py`). Centroid shift biases shapes — this matters for UNIONS shear. **Does the fix actually remove the bias?** Look for the before/after evidence; if it's not in the PR, that's the key thing to ask for (or to *run* on candide).
- **`r50` as the galaxy-size guess** (changed from `T`). Reasonable? Convergence/robustness implications for the fitter guesses.
- **New `src/shapepipe/testing/simulate.py` + `scripts/python/fitting.py`** — what do they test/simulate, and do they give a cheap way to validate the v2.0 fitter independent of a full pipeline run?

### Run it on candide (the repo is right here)

`project_dir` is the live shapepipe checkout. Don't review the diff cold — exercise it:

- `gh pr checkout 740` (or fetch `martin/ngmix_v2.0`) into a worktree so `develop`/`docs-rework` stays clean.
- Run the **ngmix module / its tests**; run `test_centroid_shift.py` and the centroid configs if feasible; check the v2.0 fitter produces sane output on a small input. Use the repo's `CLAUDE.md` for build/run conventions and the `uv` lockfile.
- Confirm **CI status** via `gh pr checks 740`.
- Capture what you ran and what you saw — that empirical check is the most valuable part of the review (the curious eye; don't report a result you didn't observe).

## Interactive — how this runs

`claude-opus`, **interactive**. Do the full autonomous prep first: read the diff (`gh pr diff 740`), read the changed `ngmix.py` / runner / Lucy's classes closely, run what you can on candide, and assemble (a) the checklist pass and (b) the science-quality read, with specific line-level notes staged as *draft* comments. **Then hold and wait for Cail to attach** — he'll talk through the science, sharpen the comments, and decide what's posted, what's pushed, and whether to approve. Posting to GitHub and any push happen *with* him, in his voice; the worker drafts, Cail directs. Do not close or `kill` until Cail signals done.

## Context

- **Repo:** `/automnt/n17data/cdaley/unions/shapepipe` on candide. Remotes: `origin` = CosmoStat/shapepipe, `martin` = martinkilbinger/shapepipe-1. `gh` is authenticated here.
- **PR #740** by Martin Kilbinger, opened 2026-06-01; head `ngmix_v2.0`, base `develop`; OPEN. Reviewer checklist is in the PR body.
- **@lbaumo = Lucy Baumont** — her ngmix classes/interface are central; her wrapper cleanup came out of her visit (see [[ngmix-update]]).
- **Prior art / related:** [[ngmix-update]] (the future-intent this realizes), [[prs-in-flight]] (the PR tracker — add #740 and its disposition), [[shapepipe]] (the Martin collaboration root). The repo `CLAUDE.md` carries build/test/run conventions — read it.
- **Authorship:** GitHub comments in Cail's voice, signed "— Claude on behalf of Cail" unless Cail asks otherwise on the call.
