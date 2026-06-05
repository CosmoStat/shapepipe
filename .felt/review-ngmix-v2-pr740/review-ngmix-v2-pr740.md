---
created-at: 2026-06-01T11:50:17.967872934+02:00
horizon: now
name: "Review + work: ngmix v2.0 (PR #740)"
outcome: |-
  INTERACTIVE — awaiting Cail. Martin engaged this morning (06-05 07:14–07:33Z), replying to 6 of the 11 findings: greenlit removing the dead CHECK_EXISTING_DIR resume (254) and making 51*51 configurable via STAMP_SIZE (766); explained the any→all zero-pixel loosening as intentional (737, "any==0 passed ~no stamps"); pointed weight-norm (949) at issue #604 (the 2023 0/1-mask → inverse-variance design discussion — context, not a direct answer to the sigma_mad-vs-v1 sub-finding); and opened the r50/T size-naming question, poking @lbaumo. I RAN Martin's explicit check: confirmed pars[4]=T=2σ² (l.911 uses Gaussian(σ=sqrt(pars[4]/2))); galaxy r50 (l.415) stores T (an area), PSF r50_PSFo (l.676) = sqrt(T/2) = σ — different scales, neither is the half-light radius 1.1774σ. Bug verified. Still unanswered (5): 293, 1140, runners, fitting.py, 1068 — and *_psfo reconvolved-PSF (1045), now the LONE clear merge-gate. report.html refreshed (Martin-disposition column + verified r50/T metric + re-sequenced steps). REPORT/ANALYSIS ONLY — no commits pushed. Held for Cail to direct: push Bucket A? r50/T transform-vs-rename? how to reply to Martin's six? Merge stays Cail's/Martin's gesture.
shuttle:
  agent: claude-opus
  enabled: true
  host: candide
  interactive: true
  kind: oneshot
  project_dir: /automnt/n17data/cdaley/unions/shapepipe
status: active
tags:
  - constitution
  - shapepipe
  - ngmix
  - review
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

## Round 2 (2026-06-05): next-steps report from Martin's responses — SUPERSEDES the merge-comment scope for THIS round

**Cail's ask (night of 2026-06-05):** read Martin's comments/responses to the standing two-part review (now on #741), do a *further code-review pass*, and produce a **report of recommended next steps** given where Martin has landed.

**Desired State (round 2):**
- Read every Martin (`martinkilbinger`) response across **#740** and **#741** — inline replies, review threads, and any new commits/force-pushes since the part-2 review (2026-06-03). Note the `# MKDEBUG` markers still in the code. (As of dispatch, Martin's engagement is light: one COMMENTED review + one inline reply at `ngmix.py:951` "Good. This probably was previous code before we changed to a random seed per tile." Verify whether more has landed.)
- **Further code-review pass** over the ngmix module + the 4 runners whose decorator contracts part-2 flagged (`sextractor_runner` et al.): for each of the 11 line-anchored findings, determine whether Martin's v2.0 now addresses it or it still stands.
- Produce a **next-steps report**: `report.html` in the fiber dir **and** a posted PR summary comment on #741 (Cail's voice / "— Claude on behalf of Cail"). For each open item: `{resolved by Martin | still open | needs a decision}`, owner (`Cail / Martin / Lucy`), and the **recommended sequence to get #740/#741 mergeable**.
- **REPORT ONLY — do NOT push new commits to the ngmix branch this round.** Pure analysis + PR comment. Keeps it collision-free while other shapepipe/pure_eb work is in flight.

**Exit:** autonomous — produce the report + PR comment, then close for Cail's review.
