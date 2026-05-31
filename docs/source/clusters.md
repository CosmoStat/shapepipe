# Running on a Cluster

ShapePipe runs the same way on every cluster: **through the container**. You
pull the slim `runtime` image once, bind-mount your clone, and run
`shapepipe_run` (or the CANFAR submission tooling) inside it — there is no
environment to install or activate on the host. This page covers the shared
pattern, then the specifics for each supported machine.

For what is *inside* the image and how it is built, see
[Container Workflow](container.md).

## The pattern

Three things hold on any cluster:

- **The container is the unit of execution.** Pull the runtime image to a SIF
  (`apptainer pull … docker://ghcr.io/cosmostat/shapepipe:<tag>-runtime`) and run
  the pipeline inside it. Nothing else is installed on the host.
- **Bind-mount your clone at the same path.** The config files reference their
  location for the input and output directories; bind-mounting the host clone at
  the *identical* path inside the container makes those paths resolve the same in
  and out of the container.
- **Keep SIFs and Apptainer's scratch off a quota-limited `$HOME`.** A pull under
  a tight home quota fails with `disk quota exceeded`; point `APPTAINER_CACHEDIR`
  (and the SIF itself) at a roomy data partition.

MPI jobs add one constraint: the container's MPI must speak the same PMIx wire
protocol as the host launcher (see [candide](#candide-slurm)).

(candide-slurm)=
## candide (SLURM)

candide uses **SLURM** (`sbatch`; the old `qsub`/PBS commands are gone). The repo
ships ready job scripts in `example/pbs/` — `candide_smp.sh` (single node,
parallelised with joblib) and `candide_mpi.sh` (multi-node, hybrid MPI). To run
the bundled single-tile example end to end:

```bash
# Keep the SIF and Apptainer's scratch off the quota-limited $HOME.
export DATA=/n17data/$USER                 # adjust to your data partition
export APPTAINER_CACHEDIR=$DATA/.apptainer

# Pull the runtime image (~850 MB).
apptainer pull "$DATA/shapepipe-runtime.sif" \
    docker://ghcr.io/cosmostat/shapepipe:develop-runtime

# Submit. SPDIR is your clone, bind-mounted at the same path inside the
# container; SP_IMAGE is the SIF. The same script serves the example and a real
# run — point the config inside it at your own pipeline.
SP_IMAGE="$DATA/shapepipe-runtime.sif" SPDIR="/path/to/shapepipe" \
    sbatch example/pbs/candide_smp.sh
```

Partitions are `comp` (2-day limit) and `compl` (5-day). Both job scripts read
`SP_IMAGE` (the SIF) and `SPDIR` (the clone) from the environment, so the only
thing that changes between the example and a real run is the config the script
points at.

For multi-node MPI, use `candide_mpi.sh`. It `module load`s a host OpenMPI in the
5.0.x series to match the PMIx that the image's OpenMPI speaks; a mismatched host
MPI (e.g. OpenMPI 4) silently degrades the job to *N* independent "rank 0 of 1"
processes instead of one *N*-rank job. The script's header comments carry the
detail.

Adapting these scripts to another SLURM cluster is mostly the `#SBATCH`
directives and the `module load` line — the `apptainer exec` invocation carries
over unchanged.

## CANFAR

CANFAR submission does not go through a batch scheduler. Instead you submit
container jobs to CANFAR's headless system with the `canfar_submit_job` console
script (backed by the `canfar` library), and watch them with `canfar_monitor` /
`canfar_monitor_log`. Pipeline steps are **bit-coded** through `-j` (the same
scheme as `scripts/sh/job_sp_canfar.bash`), the PSF model is chosen with
`-p psfex|mccd`, and `-V` selects the image version:

```bash
# Submit pipeline step(s) for the configured tiles (bit-coded -j).
canfar_submit_job -j 1 -p psfex -V 1.1

# Monitor sessions/jobs and stream logs.
canfar_monitor
canfar_monitor_log
```

The full production run — input preparation, the per-step `-j` table, and
post-processing — is documented in the
[CANFAR production walkthrough](pipeline_canfar.md).

```{note}
The CANFAR production submission scripts (`scripts/sh/job_sp_canfar*.bash`) still
run under the pre-container environment and are slated for the same
container-first cleanup the candide scripts received. Treat the walkthrough as
the current-but-evolving production procedure.
```

## ccin2p3

ccin2p3 is not yet containerised. The `example/pbs/cc_{smp,mpi}.sh` scripts target
the pre-container environment; a container-first path mirroring candide is future
work.
