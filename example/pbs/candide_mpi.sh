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

# Full path to environment
export SPENV="$HOME/.conda/envs/shapepipe"

# Full path to example config file and input data
export SPDIR="$HOME/shapepipe"

# Load modules
module load intelpython/3
module load openmpi/4.0.5

# Activate conda environment
source activate $SPENV

# Run ShapePipe using full paths to executables
$SPENV/bin/mpiexec --map-by node $SPENV/bin/shapepipe_run -c $SPDIR/example/config_mpi.ini

# Return exit code
exit 0
