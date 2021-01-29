#!/bin/bash

##########################
# SMP Script for CANDIDE #
##########################

# Receive email when job finishes or aborts
#PBS -M martin.kilbinger@cea.fr
#PBS -m ea

# Set a name for the job
#PBS -N sp_mb_j

# Join output and errors in one file
#PBS -j oe

# Set maximum computing time (e.g. 5min)
#PBS -l walltime=01:00:00

# Request number of cores
#PBS -l nodes=n01:ppn=4

# Full path to environment
export SPENV="$HOME/.conda/envs/shapepipe"
export SPDIR="$HOME/astro/repositories/github/shapepipe"

# Activate conda environment
#module load intelpython/3-2020.1
source activate $SPENV

echo "pwd=`pwd`"

survey="ps3pi_cfis"

export SP_RUN=.
export SP_CONFIG=/home/mkilbing/astro/repositories/github/shapepipe/example/$survey

echo "job=$job"

# Run ShapePipe using full paths to executable
$SPDIR/scripts/sh/job_sp_multiband_psf.sh 259.286 -s $survey -j $job

res=$?

exit $res
