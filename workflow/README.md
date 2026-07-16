# ShapePipe Snakemake orchestration

Snakemake workflow that orchestrates real-data ShapePipe runs, replacing the
`curl_canfar_local.sh → run_job_sp_canfar_v2.0.bash → job_sp_canfar_v2.0.bash`
bash layers (and per-site sbatch reimplementations). **Module code is untouched**:
rules call `shapepipe_run -c <config>` on the existing config chains. Design and
rationale: [CosmoStat/shapepipe#848](https://github.com/CosmoStat/shapepipe/issues/848)
(the living PRD).

## Quick start (nibi)

**One-allocation mode.** `sbatch` one node, then run Snakemake with the *local*
scheduler inside it: it bin-packs the DAG into the allocation (no per-job
sbatch). Each job's shell runs `shapepipe_run` **inside the container** (the
profile's apptainer software-deployment — you never type `apptainer`).

```bash
# One-time: a snakemake env on a SHARED filesystem (/project — NOT /tmp, which is
# node-local; the executor re-invokes this python inside jobs).
uv venv /project/def-mjhudson/cdaley/snakemake-env --python 3.12
source /project/def-mjhudson/cdaley/snakemake-env/bin/activate
uv pip install 'snakemake>=9,<10'

# Edit workflow/config.yaml: tile_list, run_dir, config_src, container, star_cats.

# The committed launcher loads apptainer/1.4.5 + the /project venv, so a fresh
# shell always has the right state. Inside your allocation:
workflow/bin/sp prepare   # prepare_tiles -> build_index.py -> prepare_exposures
workflow/bin/sp run       # the main compute DAG
```

`sp` subcommands: `prepare` (chains all three prepare steps), `prepare_tiles`,
`index`, `prepare_exposures`, `run`, `rerun <target>` (forced recompute of
protected products), `cancel <run>` (scancel sweep before `--unlock`),
`clean-exposures [--apply]`, `report`.

## Execution is three static invocations

The exposure job set is data-derived from the tiles' `find_exposures` output, so
it cannot be scheduled in the same static DAG that produces it. Hence:

1. `snakemake prepare_tiles` — per-tile static DAG (Git_vos → Uz → Fe),
   `keep-going` so tile failures are independent.
2. `build_index.py` — a **plain script**, not a DAG node. It iterates the
   declared tile list, checks each tile's Fe output at its deterministic path (no
   glob), writes `run_index.sqlite` + `missing.json`, and exits nonzero only if
   the missing *fraction* exceeds a threshold. (There is no `--allow-missing`
   Snakemake flag; that mechanism does not exist.)
3. `snakemake prepare_exposures` — per-exposure static DAG (Gie_vos + per-unit
   star cats), parseable now the index exists.
4. `snakemake` (main DAG) — the tile/exposure compute chain.

`sp prepare` chains 1–3 so the UX stays one command.

## Layout

```
workflow/
  Snakefile              parse-time index load; global container:; onsuccess/onerror report hooks
  config.yaml            the run: tile list, paths, container, chunk count
  bin/sp                 committed launcher (module load + /project venv + subcommands)
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
    clean_exposures.py   index-driven exposure-store reclaim (replaces temp())
profiles/nibi/config.yaml  local scheduler; apptainer SDM; resources; keep-going
```

## How it works

- **Completeness = the count floor.** Each stage declares its deterministic
  `output/run_sp_<prefix>/` directory as `directory()` output (the per-unit
  config copy forces `RUN_DATETIME = False`, so the path is known at DAG time).
  After the run, `sp_rule.py` counts products per mandatory runner against
  `completeness.py`'s floor and exits nonzero below it; per-CCD attrition between
  floor and `expect` is tolerated (the directory simply lacks that product).
  This *is* the failure policy — no 3-class taxonomy, no error-signature
  whitelist. `--keep-going` isolates a failure to its own cone.
- **Per-unit isolation is by work-dir content, not `-e/--exclusive`.** Every
  tile/exposure runs its own `shapepipe_run` in `tiles/<ID>/` or `exp/<base>/`,
  with a `cfis` config symlink, `star_cat_{exp,tiles}` symlinks, a
  `tile_numbers.txt` (tiles) or fabricated `exp_numbers-000-000.txt` pseudo-Fe
  (exposures), and a config copy setting `RUN_DATETIME=False` and
  `SMP_BATCH_SIZE={threads}`. `NUMBER_LIST=-<ID>` is injected **only** for
  stages whose numbering scheme is the unit ID (tile-scheme stages + exp split),
  never for get_images / exp_mask / exp_psf.
- **No `temp()` on shared exposure dirs.** Exposures overlap tiles, so temp()
  would cascade destructive reruns when a tile is appended. The store is
  persistent; `sp clean-exposures` reclaims space only for exposures all of
  whose consuming tiles have a final_cat.
- **The index is parse-time data, never a rule input.** Appending tiles changes
  which jobs exist without invalidating completed work. Star cats are re-keyed
  per unit for the same reason. `final_cat` is `protected()`.
- **Failure is a report, not a gate.** `run_report.py` disk-scans the trees
  against the count table and enumerates shortfalls (whole-unit absence vs
  per-CCD attrition). It is a standalone script — a DAG report node would be
  poisoned by the very failures it must enumerate — emitted automatically by the
  `onsuccess`/`onerror` hooks and runnable any time via `sp report`.
