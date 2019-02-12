#!/bin/bash

##########################
# SMP Script for ccin2p3 #
##########################

# Receive email when job finishes or aborts
#$ -M <name>@cea.fr
#$ -m ea
# Set a name for the job
#$ -N shapepipe_smp
# Set a group for the job
#$ -P P_euclid
# Join output and errors in one file
#$ -j y
# Request muliprocessing resources
#$ -l os=cl7

# Activate conda environment
ccenv anaconda
source activate $HOME/.conda/envs/shapepipe

# Run ShapePipe
cd $HOME/ShapePipe
$HOME/.conda/envs/shapepipe/bin/python shapepipe_run.py -c example/pbs/config_smp.ini

# Return exit code
exit 0
