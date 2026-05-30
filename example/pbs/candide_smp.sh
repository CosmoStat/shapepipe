#!/bin/bash

##########################
# SMP Script for CANDIDE #
##########################

# Receive email when job finishes or aborts
#PBS -M <name>@cea.fr
#PBS -m ea
# Set a name for the job
#PBS -N shapepipe_smp
# Join output and errors in one file
#PBS -j oe
# Set maximum computing time (e.g. 5min)
#PBS -l walltime=00:05:00
# Request number of cores
#PBS -l nodes=4

# Path to the local ShapePipe clone (holds the example configs and data)
export SPDIR="${SPDIR:-$HOME/shapepipe}"

# Path to the ShapePipe runtime image. Pull it once with:
#   apptainer pull "$SP_IMAGE" docker://ghcr.io/cosmostat/shapepipe:develop-runtime
export SP_IMAGE="${SP_IMAGE:-$HOME/shapepipe_develop-runtime.sif}"

# Run ShapePipe through the container -- no Python environment to activate. The
# clone is bind-mounted at the same path so that $SPDIR resolves identically
# inside the container, where the config references it for the input and output
# directories.
apptainer exec \
    --bind "$SPDIR:$SPDIR" \
    --env SPDIR="$SPDIR" \
    "$SP_IMAGE" \
    shapepipe_run -c "$SPDIR/example/pbs/config_smp.ini"

# Propagate the pipeline's exit code to the batch system
exit $?
