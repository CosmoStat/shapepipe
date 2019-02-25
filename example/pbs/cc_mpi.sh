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
# Request muliprocessing resources
#$ -l os=cl7
# Request number of cores
#$ -pe openmpi 8

# Activate conda environment
ccenv anaconda
ccenv openmpi
source activate $HOME/.conda/envs/shapepipe

# Run ShapePipe
cd $HOME/ShapePipe
/usr/local/openmpi/2.1.1/bin/mpiexec -n $NSLOTS $HOME/.conda/envs/shapepipe/bin/python shapepipe_run.py -c example/pbs/config_mpi.ini

# Return exit code
exit 0
