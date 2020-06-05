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
# Request number of cores (e.g. 40)
#PBS -l nodes=4:ppn=10
# Allocate total number of cores to variable NSLOTS
NSLOTS=`cat $PBS_NODEFILE | wc -l`

# Load moudules and activate conda environment
module load intelpython/3
module load openmpi/4.0.2
source activate $HOME/.conda/envs/shapepipe

# Run ShapePipe using full paths to executables
cd $HOME/shapepipe
/softs/openmpi/4.0.2-torque-CentOS7/bin/mpiexec -n $NSLOTS $HOME/.conda/envs/shapepipe/bin/shapepipe_run -c example/pbs/config_mpi.ini

# Return exit code
exit 0
