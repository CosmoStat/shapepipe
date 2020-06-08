#!/bin/bash

##########################
# MPI Script for ccin2p3 #
##########################

# Receive email when job finishes or aborts
#$ -M <name>@cea.fr
#$ -m bea
# Set a name for the job
#$ -N shapepipe_mpi
# Set a group for the job
#$ -P P_euclid_sci
# Join output and errors in one file
#$ -j y
# Set maximum computing time (e.g. 5min)
#$ -l h_cpu=00:05:00
# Request muliprocessing resources
#$ -l os=cl7
# Request number of cores
#$ -pe openmpi 4

# Full path to environment
export SPENV="$HOME/.conda/envs/shapepipe"
export SPDIR="$HOME/shapepipe"

# Activate conda environment
ccenv anaconda
ccenv openmpi
source activate $SPENV

# Run ShapePipe
$SPENV/bin/mpiexec -n $NSLOTS $SPENV/bin/shapepipe_run -c $SPDIR/example/pbs/config_mpi.ini

# Return exit code
exit 0
