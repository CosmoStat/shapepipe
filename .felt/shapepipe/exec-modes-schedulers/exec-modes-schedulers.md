---
name: 'ShapePipe execution modes (smp/mpi) and schedulers (PBS/SLURM): what the repo''s tooling shows'
tags:
    - shapepipe
    - mpi
    - reference
created-at: 2026-05-31T16:51:46.221097637+02:00
outcome: 'By the repo''s lights SMP is the exercised path (55/56 example configs; every canfar/candide job script is SMP-only via N_SMP, SLURM+conda); MPI is the 2019 mode, set in 1 config, and its code/config drifted out of sync (module_config_sec bug dates to #415 by git history). PBS is dead (2019 example scripts only); SLURM is current everywhere. CAVEAT: this is what the repo shows, not how ShapePipe was actually run — canfar carried most processing and is invisible from here, so MPI usage history is unknown.'
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
  example config uses it. The `worker()` call in `mpi_run.py` has been out of sync
  since PR #415 (Jan 2025) — `worker()` gained a `module_config_sec` param and
  `mpi_run.py` wasn't updated, so it passes 7 args where 8 are required. On candide
  it couldn't even wire up (PMIx mismatch, see [[shapepipe/mpi-hybrid]]), so the
  code bug couldn't surface here. Whether MPI was run elsewhere (canfar especially,
  which we can't see) is unknown — what's clear is the repo's tooling is all SMP.

**SMP and MPI are the same computation behind two dispatchers.** Both call the
identical `WorkerHandler.worker()` with the identical 8 args (`job_handler._distribute_smp_jobs`
vs `mpi_run.submit_mpi_jobs`). The MPI path's only inter-rank traffic is `bcast`
of setup objects, one `scatter` of the independent job-list, and one `gather` of
result dicts — `worker_handler.py` (the actual work) has zero MPI in it. No
`Send`/`Recv`/`Allreduce`/`Barrier` during compute. That's the signature of an
**embarrassingly parallel** workload: MPI provides no computational capability
that SMP-on-a-node-plus-a-scheduler lacks — it's a job-distribution convenience
(one `mpirun` spanning nodes vs. the submission layer fanning out per-node jobs).
This is what grounds the "is MPI worth keeping?" question to Martin — observed
from the comm pattern, not inferred from usage.

Note `MODE` is overloaded across config sections — `CLASSIC`, `MULTI-EPOCH`,
`FIT_VALIDATION`, `VALIDATION` are *module* modes (PSF / ngmix), not `[EXECUTION]`
modes. Only `smp`/`mpi` live under `[EXECUTION]`.

## Axis 2 — scheduler (the batch wrapper, outside ShapePipe)

- **PBS** (`#PBS` / `qsub`) — the 2019 `example/pbs/` scripts. **Dead** on candide
  (migrated to SLURM). All `#PBS` directives removed on the #737 branch.
- **SLURM** (`#SBATCH` / `sbatch`) — **current everywhere**. canfar since ~2020,
  candide since 2024.

## What the dates and tooling show

The maintained submission tooling is SMP-only and SLURM-based: `scripts/sh/run_scratch_local.sh`
(2024-11, *"submit jobs on candide"*) → `init_run_exclusive_canfar.sh` → `job_sp_canfar.bash`,
all `sbatch`, all **SMP** via `N_SMP` ("SMP mode only" in their help), and still **conda**
(`CONDA_PREFIX=$HOME/.conda/envs/shapepipe`), *not* the container. The `example/pbs/candide_{smp,mpi}.sh`
scripts are 2019 **teaching examples** (untouched until the #737 branch).

This is evidence about the tooling, not a claim about run history. It's suggestive — the
SMP tooling is what's been maintained, the MPI mode and its example config drifted untouched —
but most processing ran on canfar, which isn't visible from this repo, so how much MPI was
actually used is a question for the people who ran it, not something the repo can answer.

## Implications

- The MPI fix is worth landing — `mpi` is a supported mode and getting it working through
  the container on candide was the point — framed as enablement/verification, not as
  unblocking some known-active workload.
- Production scripts (SMP + SLURM + conda) are untouched by #737 and out of scope; they're
  also **not yet containerized** — a future gap to name.
- **Decision deferred to Martin (asked in #737):** is MPI worth getting working /
  maintaining on candide at all, or should candide just use SMP (which works through
  the container — `candide_smp.sh`)? Given SMP and MPI are the same computation, MPI
  earns its keep only as an ergonomic convenience. We do *not* retire it unilaterally —
  it's a documented public mode; #737 leaves it in working order and Martin makes the
  call. If kept, add a CI smoke so it can't silently rot again; if dropped, removal is
  clean and contained (`mpi_run.py`, `run_mpi`, the `import_mpi` branches, `mpi4py`,
  `candide_mpi.sh`).
