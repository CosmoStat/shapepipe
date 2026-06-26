---
id: 01KTCHWZWYYTDSYHAP4FPE7V3V
name: PRs in flight after v2 merge
tags:
    - shapepipe
    - pr
created-at: 2026-04-27T11:26:49.300097608+02:00
outcome: 'Post-v2 + post-propagation: infra stream now landed (#718 setuptools, #719 uv-lockfile, #728 dependabot+SHA-pin), supply-chain hygiene done (20 → 0 alerts). Issue #712 empirically verified resolved against current `:develop` (all 11 packages in Martin''s May 18 list import in both read-only and writable sandbox modes); comment posted, awaiting Martin reply before closing. Science PRs still open: #714 develop-bugs (closes #709 + #711 only — #712 closes separately), #708 testing-scaffold (mine); #725 centroid shift (Axel), several older Martin PRs (#704 #703 #699 #660 #650 #636), #670 lbaumo file_io. Next thread: merge #714.'
insights:
    714-already-redundant:
        claim: 'Surprise from rebasing #714: its Dockerfile commit (cf304f8f, adding astroquery/numba/fitsio + setuptools<81 pin) was *already* redundant on current develop — the v2 merge silently put astroquery/numba/fitsio into pyproject and the v2 Dockerfile installs them via ''pip install -e ".[fitsio]"'' at the end. setuptools<81 went away via #718. So ''rebase to drop the obsolete commit'' wasn''t waiting on #719 — it was already obsolete the moment v2 merged. Worth checking sooner next time before assuming a fix is still load-bearing.'
    xfail-mostly-fixable:
        claim: 'Most #708 xfails are about to be resolved: canfar_monitor IndentationError (4 xfails) and summary_run -h (1 xfail) are fixed in #714; astroquery/numba/fitsio import xfails (5 modules) resolve in #719 because uv sync installs them from pyproject. Only stile/treecorr corr2 (4 modules) is a separate issue requiring stile removal or upstream patch.'
    dependabot-policy:
        claim: 'shapepipe now ships `.github/dependabot.yml` (#728) with 14-day cooldown, monthly grouped lockfile PRs, github-actions ecosystem opted in, and SHA-pinned actions across all four workflows. Reasoning lives in the file itself + the #728 PR body. Companion fiber [[shapepipe/sqlitedict-pickle-smell]] tracks the single dismissed alert.'
    712-empirically-resolved:
        claim: 'Issue #712 is empirically resolved against current `ghcr.io/cosmostat/shapepipe:develop` (dev target, post-#728). Both the original packages (astroquery, numba, fitsio) and Martin''s May 18 follow-up list (scipy, joblib, importlib_metadata, tqdm, LSSTDESC.Coord, pyyaml, astropy_iers_data, pyerfa) import cleanly in both read-only and writable sandbox modes, as do the three originally-flagged runner modules. Pyproject confirms astroquery/numba/joblib/tqdm are core deps; the rest are transitives of astropy/mccd/modopt/galsim; fitsio is gated in both runtime (`--extra jupyter --extra fitsio`) and dev (`--extra dev`) targets. Comment posted; awaiting Martin reply before closing. Likely root cause of the May 18 report: cached/older image.'
decisions:
    setuptools-pin:
        label: drop setuptools<81 pin
        default: merged
        options:
            merged:
                label: 'Already merged as #718 (c9e71df8) — small one-liner, agreed in transcript'
---

Snapshot of CosmoStat/shapepipe PR state, maintained as a living index.

## Open — infra

(All infra PRs landed. The dependabot stream is resolved; supply-chain
posture set; SHA-pins in place. See [[shapepipe/sqlitedict-pickle-smell]]
for the one open security-fiber.)

## Open — issues (mine)

| # | What | Status |
|---|---|---|
| #712 | Dockerfile missing runtime deps | Empirically resolved against current `:develop` ([comment](https://github.com/CosmoStat/shapepipe/issues/712#issuecomment-4562085977)). Both original list (astroquery/numba/fitsio) and Martin's May 18 follow-up (scipy/joblib/importlib_metadata/tqdm/LSSTDESC.Coord/pyyaml/astropy_iers_data/pyerfa) import cleanly in read-only + writable sandbox modes. Awaiting Martin reply before closing. |
| #711 | summary_run -h crashes | Fixed by #714 (auto-closes on merge) |
| #709 | canfar_monitor IndentationError | Fixed by #714 (auto-closes on merge) |

## Open — mine (science / fixes)

| # | Branch | What | Status |
|---|---|---|---|
| #731 | `chore/smoke-test-read-only` | smoke-test in read-only mode | Open. Adds `shapepipe_run_example` wrapper; CI now runs the entry-point smoke under `docker --read-only --tmpfs /tmp:rw`. See [[shapepipe/smoke-test-read-only]]. |
| #714 | `fix/develop-bugs` | small develop bugs (#709, #711) | Open. Originally a multi-bug fix; the Dockerfile portion got absorbed into #719. Worth checking what's still load-bearing here vs already-fixed-upstream. |
| #708 | `chore/testing-scaffold` | Tier 0–2 test scaffolding | Open. Some xfails should have flipped to xpass after the v2 + uv-lockfile work; needs a rebase + xfail-list audit. |

## Open — others' PRs awaiting attention

| # | Author | What |
|---|---|---|
| #741 | martinkilbinger / lbaumo | **Ngmix v2.0** — upstream ngmix 2.4.0 + Lucy's new classes. Canonical PR (CI mirror; fork PR #740 closed by Martin). OPEN, mergeable, CI green. Two-part review + next-steps triage delivered; 11 findings open (5 cut-and-dry, 5 decisions, 1 resume), 2 are merge-gates (weight-norm, `*_psfo`). See [[review-ngmix-v2-pr740]]. |
| #725 | aguinot | Fix centroid shift (overlaps #741 centroid work — cross-ref/supersede?) |
| #704 | martinkilbinger | Contributors |
| #703 | martinkilbinger | V1.3.x |
| #699 | martinkilbinger | Coverage mask |
| #670 | lbaumo | file_io handles sextractor header |
| #660 | martinkilbinger | Existing output directory |
| #650 | martinkilbinger | Third-party catalogue for tile objects |
| #636 | martinkilbinger | Rho statistics: flexible training/test split |

## Recently closed

- **#728** `chore/dependabot-config` — dependabot.yml + SHA-pin all actions. Merged 2026-05-28.
- **#727, #726, #724, #722, #721, #720** — dependabot security bumps for idna/urllib3/gitpython/mistune/jupyter-server/jupyterlab. All squash-merged 2026-05-28 (see [[shapepipe/dependabot-pr-triage]]).
- **#719** `chore/uv-lockfile` — merged 2026-05-05 (Martin).
- **#718** `chore/drop-setuptools-pin` — merged.
- **v2.0 PR** — merged. Source of the skaha/conda situation that #719 unwound.

## Connections

- [[shapepipe]] — root
- [[shapepipe/docker-uv-revert]] — drove #719
- [[shapepipe/dependabot-pr-triage]] — drove the 6 security-bump merges (closed)
- [[shapepipe/sqlitedict-pickle-smell]] — future-work fiber for the one dismissed alert
