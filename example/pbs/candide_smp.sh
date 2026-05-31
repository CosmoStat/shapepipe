#!/bin/bash
#
# SMP (single-node) Apptainer job for candide (SLURM).
#
# SMP mode parallelises with joblib inside a single process across the allocated
# cores -- no host MPI is involved. Use this for single-node runs; use
# candide_mpi.sh to span multiple nodes.
#
# Submit with:  sbatch candide_smp.sh

#SBATCH --job-name=shapepipe_smp
#SBATCH --partition=comp
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=00:05:00
#SBATCH --output=%x-%j.log
## #SBATCH --mail-type=END,FAIL
## #SBATCH --mail-user=<name>@cea.fr

# Path to the local ShapePipe clone (holds the example configs and data).
export SPDIR="${SPDIR:-$HOME/shapepipe}"

# Path to the ShapePipe runtime image. Pull it once with:
#   apptainer pull "$SP_IMAGE" docker://ghcr.io/cosmostat/shapepipe:develop-runtime
export SP_IMAGE="${SP_IMAGE:-$HOME/shapepipe_develop-runtime.sif}"

# Run ShapePipe through the container -- no Python environment to activate. The
# clone is bind-mounted at the same path so that $SPDIR resolves identically
# inside the container, where the config references it for input / output
# directories. Keep SMP_BATCH_SIZE in config_smp.ini aligned with
# --cpus-per-task above.
apptainer exec \
    --bind "$SPDIR:$SPDIR" \
    --env SPDIR="$SPDIR" \
    "$SP_IMAGE" \
    shapepipe_run -c "$SPDIR/example/pbs/config_smp.ini"

# Propagate the pipeline's exit code to the batch system.
exit $?
