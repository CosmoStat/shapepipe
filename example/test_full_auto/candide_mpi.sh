#!/bin/bash

##########################
# MPI Script for Candide #
##########################

# Receive email when job finishes or aborts
#PBS -M axel.guinot@cea.fr
#PBS -m ea
# Set a name for the job
#PBS -N shapepipe_single_exp
# Join output and errors in one file
#PBS -j oe
# Request number of cores
#PBS -l nodes=n01:ppn=5+n02:ppn=46+n08:ppn=46
# #PBS -l walltime=00:05:00
NSLOTS=`cat $PBS_NODEFILE | wc -l`

echo $NSLOTS

# Activate conda environment
module load intelpython/3
module load openmpi/4.0.1
source activate $HOME/.conda/envs/shapepipe
export OMP_NUM_THREADS=1

# Run ShapePipe
cd $HOME/ShapePipe
/softs/openmpi/4.0.1-torque-CentOS7/bin/mpiexec -n $NSLOTS $HOME/.conda/envs/shapepipe/bin/python shapepipe_run.py -c example/test_full_auto/config_c.ini

# Return exit code
exit 0
