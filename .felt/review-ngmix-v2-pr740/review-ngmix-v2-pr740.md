---
id: 01KTCHWZX873VG28AA663A3ZQE
name: 'Review + work: ngmix v2.0 (PR #740)'
status: closed
tags:
    - constitution
    - shapepipe
    - ngmix
    - review
created-at: 2026-06-01T11:50:17.967872934+02:00
closed-at: 2026-06-05T23:59:38.502Z
outcome: 'INTERACTIVE — live with Cail. Bucket-A cleanups PUSHED to ngmix_v2.0/#741 (dd4f656a..bd60dc8e: 293/1140/254-resume/766 + fitting.py deleted); CI running. Dev image built (/n17data/cdaley/containers/shapepipe-dev) + ngmix 2.4.0 installed → full suite + example pipeline now runnable. Ran a 6+2-agent WORKFLOW on the two hard problems → two decision-ready reports in fiber: weights-report.md, size-report.md, deep-dive-report.html (sent to Cail). WEIGHTS (#604+949): found TWO coupled regressions in prepare_ngmix_weights (ngmix.py:871) — R1 noise estimator regressed from object-free windowed get_noise (now DEAD at :826) to flux-contaminated whole-stamp sigma_mad; R2 lost the v1 binarization → double-counts a real inverse-variance map. EMPIRICALLY CONFIRMED (truth ivar 1e6 → recovered 8.8e11, ratio≈1/noise²). Clean red→green test target via make_data. Rec: SPLIT — minimal v1-restore (reinstate binarization + get_noise) + test in #741; SExtractor BACKGROUND_RMS baseline as a SEPARATE PR (Codex shuttle, ~4-8h, closes #604). SIZE (r50/T): galaxy r50=pars[4]=T (area), PSF r50psf=σ; neither is 1.1774σ; all 5 r50* cols new in v2.0. UNIONS-3500 WL I (arXiv:2605.13549) reports half-light radius r_h as PRIMARY. Rec: TRANSFORM at source (honest r50 in ngmix) + cs_util converter web; bonus find — sp_validation galaxy.py:T_to_fwhm is dimensionally wrong (feeds the scale-dependent PSF-leakage fit). Both work-streams now DISPATCHED as shuttles (Cail''s go): [[ngmix-weights-ivar]] (Codex — minimal #741 regression fix w/ red→green test + the SExtractor BACKGROUND_RMS baseline as a separate PR closing #604) and [[ngmix-size-columns]] (claude-opus — honest r50 at the ngmix source + cs_util converter web + the sp_validation T_to_fwhm leakage fix). Workers prepare branches/commits/PR-descriptions + reports; merge/PR-creation stays Cail''s. STILL OPEN for the eventual #741 reply (not in these streams): *_psfo reconvolved-PSF merge-gate (1045) + the runner-decorator contract (testable via the example pipeline). PLAN: let the two shuttles finish → reassess → write the Martin reply. No reply drafted yet.'
horizon: now
shuttle:
    enabled: true
    kind: oneshot
    interactive: true
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

## Round 2 (2026-06-05): next-steps report from Martin's responses — SUPERSEDES the merge-comment scope for THIS round

**Cail's ask (night of 2026-06-05):** read Martin's comments/responses to the standing two-part review (now on #741), do a *further code-review pass*, and produce a **report of recommended next steps** given where Martin has landed.

**Desired State (round 2):**
- Read every Martin (`martinkilbinger`) response across **#740** and **#741** — inline replies, review threads, and any new commits/force-pushes since the part-2 review (2026-06-03). Note the `# MKDEBUG` markers still in the code. (As of dispatch, Martin's engagement is light: one COMMENTED review + one inline reply at `ngmix.py:951` "Good. This probably was previous code before we changed to a random seed per tile." Verify whether more has landed.)
- **Further code-review pass** over the ngmix module + the 4 runners whose decorator contracts part-2 flagged (`sextractor_runner` et al.): for each of the 11 line-anchored findings, determine whether Martin's v2.0 now addresses it or it still stands.
- Produce a **next-steps report**: `report.html` in the fiber dir **and** a posted PR summary comment on #741 (Cail's voice / "— Claude on behalf of Cail"). For each open item: `{resolved by Martin | still open | needs a decision}`, owner (`Cail / Martin / Lucy`), and the **recommended sequence to get #740/#741 mergeable**.
- **REPORT ONLY — do NOT push new commits to the ngmix branch this round.** Pure analysis + PR comment. Keeps it collision-free while other shapepipe/pure_eb work is in flight.

**Exit:** autonomous — produce the report + PR comment, then close for Cail's review.
