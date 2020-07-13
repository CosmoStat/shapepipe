#!/bin/bash

if [ -z $1 ]; then
	echo "Usage: [qsub -N <jobname>] job_PBS_SMP.sh TILE_ID"
	exit 1
fi

# Job name
#$ -N ShapePipe

# Project name
#$ -P P_euclid_ext

# Environment
#$ -pe multicores 8

#$ -l os=cl7

#$ -j y

# Mail options
#$ -M martin.kilbinger@cea.fr
#$ -m eas

echo -n "pwd = "
pwd

source ~/miniconda3/bin/activate shapepipe

# Run ShapePipe
cmd="bash $HOME/astro/repositories/github/shapepipe/scripts/sh/canfar_sp.bash $1"
echo "Running command $cmd"
$cmd

echo "Job done"

# Return exit code
exit 0
