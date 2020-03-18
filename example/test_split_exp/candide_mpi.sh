#!/bin/bash

##########################
# MPI Script for Candide #
##########################

# Receive email when job finishes or aborts
#PBS -M axel.guinot@cea.fr
#PBS -m ea
# Set a name for the job
#PBS -N shapepipe_mpi_split_exp
# Join output and errors in one file
#PBS -j oe
# Request number of cores
## -l nodes=n02:ppn=48+n03:ppn=48
#PBS -l nodes=n02:ppn=2
#PBS -l walltime=02:00:00
NSLOTS=`cat $PBS_NODEFILE | wc -l`

# Activate conda environment
module load intelpython/3
module load openmpi/4.0.0
source activate $HOME/.conda/envs/shapepipe

# Run ShapePipe
cd $HOME/ShapePipe
/softs/openmpi/4.0.0-torque-CentOS7/bin/mpiexec -n $NSLOTS $HOME/.conda/envs/shapepipe/bin/python /home/guinot/pipeline/ShapePipe/shapepipe_run.py -c /home/guinot/pipeline/ShapePipe/example/test_split_exp/config.ini

# Return exit code
exit 0
