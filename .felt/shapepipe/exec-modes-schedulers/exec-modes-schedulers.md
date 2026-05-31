---
name: 'ShapePipe execution modes (smp/mpi) and schedulers (PBS/SLURM): what''s used vs legacy'
tags:
    - shapepipe
    - mpi
    - reference
created-at: 2026-05-31T16:51:46.221097637+02:00
outcome: 'SMP is the production workhorse (55/56 example configs; all canfar/candide scripts via N_SMP, SLURM+conda); MPI is 2019 legacy used by 1 config and broken since #415. PBS is dead (2019 example scripts only); SLURM is current everywhere. The MPI module_config_sec bug survived 16mo because nobody runs MPI.'
---

Two orthogonal axes that are easy to conflate when reasoning about how ShapePipe
runs on a cluster. This fiber pins down what each is, when it entered, and what's
actually used today vs. legacy — the context for [[shapepipe/mpi-hybrid]].

## Axis 1 — execution mode (`[EXECUTION] MODE`, inside ShapePipe)

Dispatched in `src/shapepipe/run.py`: `mode = config["EXECUTION"]["MODE"].lower()`,
then `run_mpi(pipe, comm)` if `mode == "mpi"` else `run_smp(pipe)`. If mpi4py isn't
importable, mode is forced to `smp`.

- **`smp`** — joblib `Parallel(n_jobs=batch_size)` across cores on **one node**
  (`job_handler._distribute_smp_jobs`). **The living path.** 55 of 56 example
  configs set `MODE = SMP`; every canfar/candide production script drives it by
  injecting `N_SMP` into the config (`SMP_BATCH_SIZE`).
- **`mpi`** — mpi4py scatter/gather across **multiple nodes** (`pipeline/mpi_run.py`,
  `submit_mpi_jobs`). 2019-era (`c6554983` "initial mpi framework"). Exactly **1**
  example config uses it. **Broken since PR #415 (Jan 2025)**: `worker()` gained a
  `module_config_sec` param and `mpi_run.py` was never updated, so it passed 7 args
  where 8 are required. Invisible for 16 months because nobody runs MPI — and on
  candide it couldn't even wire up (PMIx mismatch, see [[shapepipe/mpi-hybrid]]),
  which masked the code bug underneath.

Note `MODE` is overloaded across config sections — `CLASSIC`, `MULTI-EPOCH`,
`FIT_VALIDATION`, `VALIDATION` are *module* modes (PSF / ngmix), not `[EXECUTION]`
modes. Only `smp`/`mpi` live under `[EXECUTION]`.

## Axis 2 — scheduler (the batch wrapper, outside ShapePipe)

- **PBS** (`#PBS` / `qsub`) — the 2019 `example/pbs/` scripts. **Dead** on candide
  (migrated to SLURM). All `#PBS` directives removed on the #737 branch.
- **SLURM** (`#SBATCH` / `sbatch`) — **current everywhere**. canfar since ~2020,
  candide since 2024.

## The story the dates tell

ShapePipe shifted from **"a few big MPI jobs under PBS on candide" (2019)** to
**"many small SMP jobs under SLURM" (2024+)**. Today's production submission path
is `scripts/sh/run_scratch_local.sh` (2024-11, *"submit jobs on candide"*) →
`init_run_exclusive_canfar.sh` → `job_sp_canfar.bash`: all `sbatch` (SLURM), all
**SMP** mode via `N_SMP`, and still **conda** (`CONDA_PREFIX=$HOME/.conda/envs/shapepipe`),
*not* the container.

The `example/pbs/candide_{smp,mpi}.sh` scripts are 2019 **teaching examples**
(untouched until #737 branch), not the production path.

## Implications

- The MPI bug fix is worth landing — `mpi` is still a supported mode and fixing it
  on candide was the goal — but it restores a *legacy* path, it doesn't unblock
  production.
- Production canfar/candide scripts (SMP + SLURM + conda) are untouched by #737 and
  out of scope; they're also **not yet containerized** — a future gap to name.
- **Open question for Martin / the team:** does anyone still need MPI multi-node
  runs at all, or has SMP-under-SLURM (many per-node jobs) fully replaced it? If MPI
  is truly dead, the honest move might be to retire it rather than maintain it.
