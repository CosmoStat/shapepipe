#!/bin/bash

##########################
# SMP Script for Candide #
##########################

# Receive email when job finishes or aborts
#PBS -M <name>@cea.fr
#PBS -m ea
# Set a name for the job
#PBS -N shapepipe_smp
# Join output and errors in one file
#PBS -j oe
# Request number of cores
#PBS -l nodes=4

# Activate conda environment
module load intelpython/3
source activate $HOME/.conda/envs/shapepipe

# Run ShapePipe
cd $HOME/ShapePipe
$HOME/.conda/envs/shapepipe/bin/python shapepipe_run.py -c example/pbs/config_smp.ini

# Return exit code
exit 0
