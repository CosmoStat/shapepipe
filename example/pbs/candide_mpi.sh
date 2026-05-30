#!/bin/bash

##########################
# MPI Script for CANDIDE #
##########################

# Receive email when job finishes or aborts
## #PBS -M <name>@cea.fr
## #PBS -m ea

# Set a name for the job
#PBS -N shapepipe_mpi

# Join output and errors in one file
#PBS -j oe

# Set maximum computing time (e.g. 5min)
#PBS -l walltime=00:05:00

# Request number of cores (e.g. 2 from 2 different machines)
#PBS -l nodes=2:ppn=2

# Path to the local ShapePipe clone (holds the example configs and data)
export SPDIR="${SPDIR:-$HOME/shapepipe}"

# Path to the ShapePipe runtime image. Pull it once with:
#   apptainer pull "$SP_IMAGE" docker://ghcr.io/cosmostat/shapepipe:develop-runtime
export SP_IMAGE="${SP_IMAGE:-$HOME/shapepipe_develop-runtime.sif}"

# Load the host MPI. ShapePipe runs as a standard "hybrid" Apptainer MPI job:
# the host mpiexec launches one container rank per slot and the in-image
# mpi4py / OpenMPI handle the communication. The image ships OpenMPI 4.1.x, so
# load a host OpenMPI in the same family for ABI compatibility.
module load openmpi

# Run ShapePipe through the container -- no Python environment to activate. The
# clone is bind-mounted at the same path so that $SPDIR resolves identically
# inside the container, where the config references it for the input and output
# directories.
mpiexec --map-by node \
    apptainer exec \
        --bind "$SPDIR:$SPDIR" \
        --env SPDIR="$SPDIR" \
        "$SP_IMAGE" \
        shapepipe_run -c "$SPDIR/example/pbs/config_mpi.ini"

# Propagate the pipeline's exit code to the batch system
exit $?
