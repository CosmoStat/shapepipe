---
name: ShapePipe hybrid MPI through the container on candide
status: active
tags:
    - shapepipe
    - mpi
    - container
    - candide
created-at: 2026-05-31T12:22:50.017370879+02:00
outcome: |-
    Two independent bugs, both fixed, verified e2e on candide. (1) LAUNCHER: container
    shipped OpenMPI 4.1.4/PMIx2 vs candide host 5.0.x/PMIx5 → hybrid MPI gave N rank-0
    singletons. Fixed by building OpenMPI 5.0.8 from source in the image (--disable-dlopen,
    bundled PMIx5/PRRTE), dropping libopenmpi-dev, keeping the mpi4py wheel (uv.lock
    untouched); SLURM-ified candide scripts; CI now publishes on every branch push.
    (2) SHAPEPIPE CODE: with ranks finally wired up, shapepipe_run under MPI hit
    "worker() missing module_runner" — a latent bug since #415 (mpi_run.py never updated
    when worker() gained module_config_sec), invisible for 16mo because MPI is the legacy
    path (SMP is production). Fixed in e5999733. Re-verified (job 780655, host-src bind
    over PR image): 4 ranks n23+n25, all 3 modules ran, real RUN_EXIT=0, 0 errors.
    REMAINING: rebuild published image with the code fix (push→CI), then Martin review +
    merge of #737.
---

## The problem

The "MPI verification gap" flagged in [[shapepipe/cleanup-rhostats-jobscripts]]:
PR #737's `candide_mpi.sh` uses the correct Apptainer **hybrid** pattern (host
`mpirun` launches one container rank per task) but couldn't be verified, and the
container/host OpenMPI versions had drifted apart.

Goal: actually run ShapePipe through the container under MPI on candide, end to
end, following [Apptainer's MPI guidance](https://apptainer.org/docs/user/main/mpi.html).

## What the data said

Empirical test on candide (image = `ghcr.io/cosmostat/shapepipe:develop-runtime`,
host `module load openmpi/5.0.8`, single node, 4 ranks):

```
mpirun -n 4 apptainer exec $SIF python -m mpi4py.bench helloworld
  → Hello, World! I am process 0 of 1 on n23.   (×4)
```

Four singletons instead of one 4-rank job. Apptainer's docs name this exactly:
*"If your containers run N rank 0 processes … the MPI stack used to launch is not
compatible with the MPI stack in the container."*

**Root cause — PMIx wire mismatch.** The hybrid model needs the container's MPI
to speak the same PMIx as the host launcher.

| | OpenMPI | PMIx |
|---|---|---|
| container (Debian bookworm `libopenmpi-dev`) | 4.1.4 | 2.x (`MCA pmix: ext3x`, `--with-pmix=.../pmix2`) |
| candide host (`openmpi/5.0.8`) | 5.0.8 | 5.x (internal) |

PMIx 2 client cannot connect to the PMIx 5 server PRRTE stands up, so each rank
initializes standalone. (`libmpi.so.40` is ABI-stable across OpenMPI 4↔5, which
is why mpi4py *imports* fine — but import isn't wire-up.)

## The fix

Build **OpenMPI 5.0.x from source** in the image (bundled PMIx 5 / PRRTE,
`--with-pmix=internal --with-prrte=internal --with-hwloc=internal
--with-libevent=internal --disable-dlopen`). The stock mpi4py wheel (from
uv.lock) dlopens `libmpi.so.40`, the soname this build provides, so it needs
**no rebuild** and `uv.lock` stays a pure SSOT. `--disable-dlopen` links MCA
components statically — it both fixes an internal-openpmix `pdl` configure
failure (wants libltdl headers otherwise) and is the right posture for a
container (no dlopen of plugin .so across the SIF/bind boundary).

Proven locally on candide before committing: a minimal proof container compiled
OpenMPI 5.0.8 + built mpi4py clean, and the `--disable-dlopen` flag was found by
iterating the configure step. Then switched to the **build-remotely / pull-
locally** loop (now in CLAUDE.md): edit Dockerfile → push → CI builds and
publishes to GHCR → `apptainer pull` on the cluster → test. Local `apptainer
build` is the wrong default — cluster quotas are tight (hit `disk quota
exceeded` on `$HOME`; keep SIFs + `APPTAINER_TMPDIR`/`CACHEDIR` on a data
partition). CI now publishes on every branch push (not just integration
branches) so any PR has a pullable, cluster-testable image before merge.

## Keeping host ↔ container MPI in sync (design)

The container seals off the host's userspace *except* MPI — to use the
interconnect + launcher you need the in-image MPI to cooperate with host
machinery you can't seal off. The contract is narrower than "same version":
what must match is the **PMIx wire protocol** and **launch mechanism**, and
PMIx is compatible *within a major version*. So the compatibility unit is the
**5.0.x series**, not the point release — hence `module load openmpi` (default)
in the job script and `OMPI_VERSION` as a Docker `ARG` (retarget = one number).

Spectrum for multi-cluster / differing-MPI futures, cheapest → most robust:
1. **Pin a series + track targets** (chosen). One image covers every PMIx-5
   cluster. Most modern HPC is here now.
2. **CI matrix → variants** from the same build-arg (`:…-ompi5`, `:…-ompi4`)
   when two targets straddle a PMIx major. One source, N artifacts.
3. **Bind model** (`--bind $MPI_DIR`): no MPI baked, host MPI mounted in —
   always matches but fragile (glibc/path/admin-bind caveats). Fallback.
4. **Wi4MPI** (a CEA tool): MPI translation layer, write-once-run-anywhere
   across MPI families. Heaviest; the escalation if 1–2 don't suffice.
5. **Preflight self-check** (complements any): run a 2-rank helloworld, detect
   the "rank 0 of 1" singleton signature, fail loudly instead of silently
   running N independent copies → wrong science. Recommended regardless; turns
   silent desync into an obvious error. Not yet implemented — candidate for
   this PR or a follow-up.

## Environment facts (candide, 2026-05)

- **Scheduler is SLURM**, not PBS — `qsub`/`qstat` are gone; partitions `comp`
  (2-day) / `compl` (5-day), idle nodes available. The `#PBS` directives in the
  candide job scripts are dead.
- **Host OpenMPI**: modules `openmpi/5.0.3`–`5.0.10`, built `-slurm-CentOS8`
  (`/softs/openmpi/5.0.8-slurm-CentOS8`). The 4.0.5 the old script loaded is gone.
- **srun launch is not viable** for OpenMPI 5 here: `srun --mpi=list` →
  none/cray_shasta/pmi2 only (no pmix). Use `mpirun` (PRRTE carries PMIx).
- **Local container builds work** via `apptainer build --fakeroot` even without
  `/etc/subuid` entries (root-mapped namespace; `allow setuid = yes`).

## Deliverables (on #737 branch `cleanup/candide-scripts-container`)

All committed (`4fc948db` MPI fix, `d31d4d26` CI), pushed, CI building. Going
onto the existing #737 PR rather than a new one — this completes the candide-
scripts work #737 started.

1. **Dockerfile** → OpenMPI 5.0.8 from source, `--disable-dlopen`; libopenmpi-dev
   dropped; mpi4py wheel kept (uv.lock untouched).
2. **candide job scripts** → SLURM (`#SBATCH`), `module load openmpi` (default),
   `mpirun -n $SLURM_NTASKS apptainer exec … shapepipe_run`.
   (`example/pbs/config_mpi.ini` already existed and is correct.)
3. **docs / CLAUDE.md** — hybrid-MPI run pattern; build-remotely/pull-locally loop.
4. **CI** — publish on every branch push so PR images are cluster-testable.
5. **ShapePipe MPI code fix** (`e5999733`) — thread `module_config_sec` through
   `run_mpi`/`submit_mpi_jobs`/`worker()`; the latent #415 bug surfaced once the
   launcher worked. Needs an image rebuild to ship.

## Empirical close (2026-05-31) — two layers

The fix turned out to have **two independent layers**. The launcher fix
(above) was necessary but not sufficient: making the ranks actually wire up
exposed a second, latent bug in ShapePipe's own MPI code.

**Layer 1 — launcher (PMIx), verified.** Pulled the PR image on candide and
ran the rank wire-up check (2 nodes, 4 tasks, `module load openmpi` → `mpirun
-n 4 apptainer exec … python -m mpi4py.bench helloworld`):

```
Hello, World! I am process 0 of 4 on n23.
Hello, World! I am process 1 of 4 on n23.
Hello, World! I am process 2 of 4 on n25.
Hello, World! I am process 3 of 4 on n25.
```

One 4-rank job spanning two nodes — the exact inverse of the pre-fix 4×
"rank 0 of 1". Image reports `Open MPI: 5.0.8`. ✓

**Layer 2 — ShapePipe MPI code, was broken, now fixed.** With the ranks wired
up, the actual `shapepipe_run` under MPI immediately hit:

```
ERROR: WorkerHandler.worker() missing 1 required positional argument: 'module_runner'
```

A latent bug since PR #415: `worker()` gained a `module_config_sec` parameter
and `pipeline/mpi_run.py:submit_mpi_jobs` was never updated, so it passed 7
args where 8 are required. Invisible for 16 months because **nobody runs MPI**
— SMP is the production path (see [[shapepipe/exec-modes-schedulers]]) and the
PMIx mismatch meant MPI never even started on candide. Fixed by threading
`module_config_sec` through `run_mpi` → `submit_mpi_jobs` → `worker()` (commit
`e5999733`), matching the SMP/serial call sites.

**Re-verified end to end** (job 780655, PR image with the working-tree `src`
bind-mounted over `/app/src` so the fix is exercised without an image rebuild):
fixed `submit_mpi_jobs` signature live in-container, 4 ranks across n23+n25,
all three modules (`python`/`serial`/`execute_example_runner`) produced output
trees, real `RUN_EXIT=0`, and `shapepipe.log` records *"A total of 0 errors
were recorded."* **Now genuinely verified.**

> Correction: an earlier close claimed the full pipeline ran clean at this
> point. It did not — that run hit the Layer-2 error and the sbatch script's
> `RUN_EXIT=0` was a hardcoded `echo`, not the real exit code. The launcher
> half was real; the pipeline half was not, until the code fix above.

**Remaining:** bake the code fix into the published image (push → CI rebuild
of `:cleanup-candide-scripts-container-runtime`), then Martin's review + merge
of #737.

(Note: the in-image `mpi4py` import looks absent under `bash -lc` because the
login shell resets PATH off the venv — a probe artifact, not real; the actual
`mpirun apptainer exec python -m mpi4py.bench` run resolves it via the image's
default PATH and wires up fine, as the helloworld output shows.)
