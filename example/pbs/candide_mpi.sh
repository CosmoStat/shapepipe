#!/bin/bash
#
# Hybrid Apptainer MPI job for candide (SLURM).
#
# ShapePipe runs as a standard Apptainer "hybrid" MPI job: the host `mpirun`
# launches one container rank per SLURM task, and the OpenMPI + mpi4py inside
# the image handle the communication. For the ranks to find one another, the
# container's OpenMPI must speak the same PMIx as the host launcher -- the
# published image ships OpenMPI 5.0.x to match candide's OpenMPI 5.0.x modules.
# (An OpenMPI 4 image silently degrades to N independent "rank 0 of 1"
# processes.)
#
# Submit with:  sbatch candide_mpi.sh

#SBATCH --job-name=shapepipe_mpi
#SBATCH --partition=comp
#SBATCH --nodes=2
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=2
#SBATCH --time=00:05:00
#SBATCH --output=%x-%j.log
## #SBATCH --mail-type=END,FAIL
## #SBATCH --mail-user=<name>@cea.fr

# Path to the local ShapePipe clone (holds the example configs and data).
export SPDIR="${SPDIR:-$HOME/shapepipe}"

# Path to the ShapePipe runtime image. Pull it once with:
#   apptainer pull "$SP_IMAGE" docker://ghcr.io/cosmostat/shapepipe:develop-runtime
export SP_IMAGE="${SP_IMAGE:-$HOME/shapepipe_develop-runtime.sif}"

# Host MPI. The image ships OpenMPI 5.0.x, and any host OpenMPI in the 5.0.x
# family is PMIx-compatible with it, so the cluster default is fine. If candide's
# default ever moves to a different major series, pin a 5.0.x here instead
# (`module load openmpi/5.0.x`) to keep the host/container PMIx match.
module load openmpi

# `mpirun` inherits the node / task layout from the SLURM allocation; -n is the
# total task count. The clone is bind-mounted at the same path so that $SPDIR
# resolves identically inside the container, where the config references it for
# the input and output directories.
mpirun -n "$SLURM_NTASKS" \
    apptainer exec \
        --bind "$SPDIR:$SPDIR" \
        --env SPDIR="$SPDIR" \
        "$SP_IMAGE" \
        shapepipe_run -c "$SPDIR/example/pbs/config_mpi.ini"

# Propagate the pipeline's exit code to the batch system.
exit $?
