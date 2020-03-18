#!/bin/bash

##########################
# SMP Script for Candide #
##########################

# Receive email when job finishes or aborts
#PBS -M axel.guinot@cea.fr
#PBS -m ea
# Set a name for the job
#PBS -N shapepipe_smp_all_tile
# Join output and errors in one file
#PBS -j oe
# Request number of cores
#PBS -l nodes=n02:ppn=10
#PBS -l walltime=100:00:00

# Activate conda environment
module load intelpython/3
source activate $HOME/.conda/envs/shapepipe

# Run ShapePipe
#cd $HOME/ShapePipe
$HOME/.conda/envs/shapepipe/bin/python /home/guinot/pipeline/ShapePipe/shapepipe_run.py -c /home/guinot/pipeline/ShapePipe/example/test_full_auto/config_c.ini

# Return exit code
exit 0
