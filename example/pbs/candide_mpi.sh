#!/bin/bash

##########################
# MPI Script for Candide #
##########################

# Receive email when job finishes or aborts
#PBS -M <name>@cea.fr
#PBS -m ea
# Set a name for the job
#PBS -N shapepipe_mpi
# Join output and errors in one file
#PBS -j oe
# Request number of cores
#PBS -l select=4:ncpus=1
NSLOTS=`cat $PBS_NODEFILE | wc -l`

# Activate conda environment
module load intelpython/3
source activate $HOME/.conda/envs/shapepipe

# Run ShapePipe
cd $HOME/ShapePipe
mpiexec -n $NSLOTS $HOME/.conda/envs/shapepipe/bin/python shapepipe_run.py -c example/pbs/config_mpi.ini

# Return exit code
exit 0
