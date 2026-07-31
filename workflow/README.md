# ShapePipe Snakemake orchestration

Snakemake workflow that orchestrates real-data ShapePipe runs. It replaces the
`curl_canfar_local.sh → run_job_sp_canfar_v2.0.bash → job_sp_canfar_v2.0.bash`
bash layers and the per-site sbatch reimplementations. **Module code is
untouched**: rules call `shapepipe_run -c <config>` on the existing config
chains. Design and rationale:
[CosmoStat/shapepipe#848](https://github.com/CosmoStat/shapepipe/issues/848)
(the living PRD).

Use `workflow/bin/sp` for everything. Bare `snakemake all` outside `sp` is
unsupported: `sp` sets the state directory, the SLURM profile, and
`SP_PHASE`, which the Snakefile needs to build the tile/exposure index at
parse time. Running snakemake directly skips all of that.

## Quick start (nibi)

```bash
# One-time: a snakemake env on a SHARED filesystem (/project — NOT /tmp, which
# is node-local; the SLURM executor re-invokes this python inside every job).
uv venv /project/def-mjhudson/cdaley/snakemake-env --python 3.12
source /project/def-mjhudson/cdaley/snakemake-env/bin/activate
uv pip install 'snakemake>=9,<10' 'snakemake-executor-plugin-slurm>=2.7,<3'

# Edit workflow/config.yaml: tile_list, run_dir, container, star_cats.

# The committed launcher loads apptainer/1.4.5 + the /project venv, so a
# fresh shell always has the right state.
workflow/bin/sp run       # bring products on disk up to date with the tile list
workflow/bin/sp report    # emit run_report.json now (mid-run is fine)
workflow/bin/sp cancel <run-name-substring>   # scancel this workflow's jobs
```

Installed and pinned versions on nibi (`/project/def-mjhudson/cdaley/snakemake-env`,
queried 2026-07-30): `snakemake==9.23.1`, `snakemake-executor-plugin-slurm==2.7.1`.
Pin range: `snakemake>=9,<10`, `snakemake-executor-plugin-slurm>=2.7,<3`. The
v8→v9 breaks matter here: `--use-singularity` became `--sdm`, executors became
plugins, and full `rerun-triggers` became the default.

Anything other than `run`, `report`, `cancel` passes straight through to
snakemake with the workflow's profile and state dir — the escape hatch for
`sp --unlock`, `sp --dag`, `sp exp_psf ...`.

## Execution: two static invocations

The exposure job set is data-derived from the tiles' `find_exposures` output,
so it cannot live in the same static DAG that produces it. `sp run` is
therefore two snakemake invocations over one Snakefile:

1. **PREPARE** — `snakemake prepare_all_tiles`: per-tile static DAG
   (`Git_vos → Uz → Fe`), `keep-going` so tile failures are independent. A
   nonzero exit here is not fatal to the run — tiles that lost their
   exposure list are dropped at the compute parse — but it is a warning:
   `SP_MISSING_THRESHOLD` (default 0.0) is the real gate.
2. **COMPUTE** — `snakemake all`: this invocation's *parse* builds the
   tile↔exposure index (`build_index.py`, imported at parse time, not a DAG
   node) and runs the full tile/exposure compute chain. The index
   accumulates across invocations, so appending tiles later changes which
   jobs exist without invalidating completed work.

`sp run` chains both so the UX is one command; both exit codes are checked
and the run fails if either phase failed.

## Execution mode: one SLURM job per rule

The profile (`profiles/nibi/config.yaml`) sets `executor: slurm`. Every rule
instance becomes its own SLURM job carrying that rule's own attempt-scaled
resources (`cpus_per_task = threads`, `mem_mb`, `runtime`) and runs inside
`apptainer exec` via the profile's software-deployment method — the workflow
never calls apptainer directly. Snakemake feeds the queue as jobs finish, so
the full campaign's job count never needs to queue at once, and multi-node
scaling is inherent. `jobs:` in the profile caps concurrent submissions at a
fraction of the cluster's per-user submit limit (queried, not invented — see
the comment in the profile file for the query and date).

The profile intentionally sets no `set-resources` / `set-threads` overrides:
those replace a rule's own values wholesale, which would kill the
attempt-scaled `mem_mb = lambda wc, attempt: ...` OOM retries and the tuned
ngmix thread count. Rules own their own resources; the profile only supplies
defaults for rules that state nothing.

`group:` labels that fuse short rules (uncompress, merges) into their chunky
neighbours (queue-latency amortization, per the PRD) are **not yet wired**:
they require labels in `workflow/rules/*.smk`, out of scope for this
profile-only pass.

## Layout

```
workflow/
  Snakefile              parse-time index load; global container:; onsuccess/onerror report hooks
  config.yaml            the run: tile list, paths, container, chunk count
  bin/sp                 committed launcher (module load + /project venv + run/report/cancel)
  rules/
    prepare.smk          tile get_images/uncompress/find_exposures + per-unit star cats
    exposure.smk         per-exposure: get_images, split, mask, psf (no temp())
    tile.smk             per-tile: exp forest, merge_headers, mask, detect, vignets, ngmix, merge, make_cat
  scripts/
    sp_rule.py           the thin per-unit wrapper (isolation furniture, config copy, log-sync, count floor)
    build_index.py       prepare-phase run_index.sqlite builder (plain script)
    build_forest.py      per-tile exposure symlink forest (group-compatible shell)
    completeness.py      the ported count-floor table (shared by sp_rule + run_report)
    run_report.py        standalone report (NOT a DAG node; run_report hooks call it)
    clean_exposure.py    ONE exposure's store + manifests -> tombstone (the clean_exposure rule)
profiles/nibi/config.yaml  SLURM executor; apptainer SDM; per-user jobs cap; keep-going
```

## How it works

- **The atom is one rule == one `shapepipe_run` on one unit.** Its single
  declared output is that unit's manifest
  (`<unit dir>/manifests/<stage>.json`), not its product files — a missing
  CCD is often legitimate, and at DR6 scale per-CCD declaration means
  millions of paths.
- **Manifests are the DAG's currency.** `completeness.py check` writes the
  manifest, and it is the only record of *why* a unit failed. That is why
  the profile sets `keep-incomplete: true` — Snakemake's default deletes a
  failed job's declared outputs, which would erase the manifest the report
  needs.
- **Completeness is a count floor, not a taxonomy.** After a run,
  `sp_rule.py` counts products per mandatory runner against
  `completeness.py`'s floor and exits nonzero below it. Per-CCD attrition
  between floor and `expect` is tolerated. No 3-class taxonomy, no
  error-signature whitelist. `--keep-going` isolates a failure to its own
  DAG cone.
- **Stores are sharded.** Every tile/exposure runs its own `shapepipe_run`
  in `tiles/<2-char prefix>/<ID>/` or `exp/<prefix>/<base>/`, with a `cfis`
  config symlink, `star_cat_{exp,tiles}` symlinks, and a config copy setting
  `RUN_DATETIME=False`. Configs are committed under `workflow/config/cfis/`
  and version with the rules that set the env vars they interpolate — there
  is no `config_src` knob, by construction.
- **The index is parse-time data, never a rule input.** Appending tiles
  changes which jobs exist without invalidating completed work. Star cats
  are re-keyed per unit for the same reason. `final_cat` is `protected()`.
- **Exposure products are not `temp()`.** Exposures overlap tiles, so
  `temp()` would cascade destructive reruns when a tile is appended later.
  Reclamation is the in-DAG `clean_exposure` rule instead: one job per
  exposure, taking every consuming tile's `tile_vignets` manifest as input
  (the campaign-wide consumer set comes from the accumulating index), which
  deletes the store *and* the exposure's manifests and leaves a `cleaned.json`
  tombstone. Deleting the manifests is what makes a late append correct: the
  appended tile finds an unbuilt chain and regenerates it. The `clean:` flag
  in `config.yaml` gates it; flipping it on later reclaims retroactively,
  since the missing tombstones schedule exactly the outstanding clean jobs.
- **Failure is a report, not a gate.** `run_report.py` disk-scans the trees
  against the count table and enumerates shortfalls (whole-unit absence vs
  per-CCD attrition). It runs standalone — a DAG report node would itself be
  poisoned by the failures it must enumerate — and fires automatically from
  the COMPUTE invocation's `onsuccess`/`onerror` hooks, or on demand via
  `sp report`.
