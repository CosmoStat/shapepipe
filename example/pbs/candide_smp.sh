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

# Activate conda environment
module load intelpython/3
source activate $HOME/.conda/envs/shapepipe

# Run ShapePipe
cd $HOME/shapepipe
$HOME/.conda/envs/shapepipe/bin/shapepipe_run -c example/pbs/config_smp.ini

# Return exit code
exit 0
