---
name: 'Review + work: ngmix v2.0 (PR #740)'
status: active
tags:
    - constitution
    - shapepipe
    - ngmix
    - review
created-at: 2026-06-01T11:50:17.967872934+02:00
outcome: 'Interactive review of PR #740 ngmix v2.0, working on our mirror PR #741. CI now FULLY GREEN (260 passed). Landed: rng-threading fix + reproducibility property test (observed red->green), scripts-import smoke test extending coverage to scripts/ (develop''s test_imports only walks the package — that''s why get_guess slipped through), and deleted the stranded v1 centroid_bias.py (old-API, imports removed get_guess, only ever run from the test_centroid_bug worktree; superseded by centroid_bias_v2.py). Centroid bug characterized: v1 harness documents m~-2.8e-2 on old ngmix 1.3.6, ~0 with fix_jac_centroid — the bug is an OLD-path phenomenon, which is why v2.0 (a full rewrite) shows no bias beyond ~+4e-4 residual. Deleted the over-scoped separate_infra_prs memory (Cail: overkill/limiting). STILL OPEN for the review write-up: PSF *_psfo doc-vs-code contradiction (now reports metacal-reconvolved PSF not original — post as question to Martin); the 4 methodology findings (weight normalization, per-type PSF collapse, zero-pixel any->all, PSF fit prior/flux) — deeper joint walkthrough with Cail NEXT; runner-decorator mention; r50/T units; harness-scope question (run_all.sh + bug/fix configs hardcode mkilbing home paths — ship only reusable v2 parts or mark local?); dangling centroid_bug/fix.yaml refs to deleted script. Nothing reviewed-posted to GitHub yet (only the #741 mirror PR body, which credits Martin).'
horizon: now
shuttle:
    enabled: true
    kind: oneshot
    interactive: true
    host: candide
    project_dir: /automnt/n17data/cdaley/unions/shapepipe
    agent: claude-opus
    session:
        id: be6fbfa6-8418-490b-bb05-d1936d5059ef
        agent: claude-opus
        dispatched_at: 2026-06-01T09:51:08.677418931Z
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
