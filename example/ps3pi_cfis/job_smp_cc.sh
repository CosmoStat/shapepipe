#!/bin/bash

##########################
# SMP Script for CANDIDE #
##########################

# Receive email when job finishes or aborts
#$ -M martin.kilbinger@cea.fr
#$ -m ea

# Set a name for the job
#$ -N sp_mb_j

# Join output and errors in one file
#$ -j y

# Set maximum computing time (e.g. 5min)
#$ -l h_cpu=05:00:00

# Request muliprocessing resources
#$ -l os=cl7

# Request number of cores
#$ -pe multicores 8

# Full path to environment
export SPENV="$HOME/miniconda3/envs/shapepipe"
export SPDIR="$HOME/astro/repositories/github/shapepipe"

# Activate conda environment
source activate $SPENV

echo "pwd=`pwd`"

survey="ps3pi_cfis"

export SP_RUN=.
export SP_CONFIG=$HOME/astro/repositories/github/shapepipe/example/$survey

# Run ShapePipe using full paths to executable
$SPDIR/scripts/sh/job_sp_multiband_psf.sh 259.286 -s $survey -j $1

res=$?

exit $res
