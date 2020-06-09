#!/bin/bash

##########################
# MPI Script for CANDIDE #
##########################

# Receive email when job finishes or aborts
#PBS -M <name>@cea.fr
#PBS -m ea
# Set a name for the job
#PBS -N shapepipe_mpi
# Join output and errors in one file
#PBS -j oe
# Set maximum computing time (e.g. 5min)
#PBS -l walltime=00:05:00
# Request number of cores (e.g. 4 from 2 different machines)
#PBS -l nodes=2:ppn=2
# Allocate total number of cores to variable NSLOTS
NSLOTS=`cat $PBS_NODEFILE | wc -l`

# Full path to environment
export SPENV="$HOME/.conda/envs/shapepipe"
export SPDIR="$HOME/shapepipe"

# Load moudules and activate conda environment
module load intelpython/3
module load openmpi/4.0.2
source activate $SPENV

# Run ShapePipe using full paths to executables
$SPENV/bin/mpiexec -n $NSLOTS $SPENV/bin/shapepipe_run -c $SPDIR/example/pbs/config_mpi.ini

# Return exit code
exit 0
