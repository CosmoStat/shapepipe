#!/bin/bash

if [ -z $1 ]; then
	echo "Usage: [qsub -N <jobname>] job_PBS_SMP.sh config.<module>.ini>"
	exit 1
fi

# Job name
#$ -N ShapePipe

# Project name
#$ -P P_euclid_sci

# Environment
#$ -pe multicores 16

#$ -l os=cl7

#$ -j y

# Mail options
#$ -M martin.kilbinger@cea.fr
#$ -m eas

echo -n "pwd = "
pwd

# Activate conda environment
# The following is from Sam's example, does not work for Martin
# (some library error)
#ccenv anaconda
#source activate $HOME/.conda/envs/shapepipe

source ~/miniconda3/bin/activate shapepipe

# Run ShapePipe
#cmd="$HOME/.conda/envs/shapepipe/bin/python $HOME/ShapePipe/shapepipe_run.py -c $1"
cmd="$HOME/miniconda3/envs/shapepipe/bin/python $HOME/ShapePipe/shapepipe_run.py -c $1"
echo "Running command $cmd"
$cmd

echo "Job done"

# Return exit code
exit 0
