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

# Full path to environment
export SPENV="$HOME/.conda/envs/shapepipe"
export SPDIR="$HOME/shapepipe"

# Activate conda environment
module load intelpython/3
source activate $SPENV

# Run ShapePipe using full paths to executables
$SPENV/bin/shapepipe_run -c $SPDIR/example/pbs/config_smp.ini

# Return exit code
exit 0
