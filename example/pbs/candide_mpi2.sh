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
#PBS -l nodes=n01:ppn=10+n02:ppn=10
# +n09:ppn=46
#PBS -l walltime=00:05:00
NSLOTS=`cat $PBS_NODEFILE | wc -l`

# Activate conda environment
module load intelpython/3
module load openmpi/4.0.0
source activate $HOME/.conda/envs/shapepipe
export OMP_NUM_THREADS=1

# Run ShapePipe
cd $HOME/ShapePipe
/softs/openmpi/4.0.0-torque-CentOS7/bin/mpiexec -n $NSLOTS $HOME/.conda/envs/shapepipe/bin/python shapepipe_run.py -c example/pbs/config_mpi.ini

# Return exit code
exit 0
