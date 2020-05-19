#!/bin/bash

##########################
# MPI Script for Candide #
##########################

# Receive email when job finishes or aborts
# #PBS -M axel.guinot@cea.fr
#PBS -m ea
# Configure job_array
#PBS -t 6501-7000%200
# Set a name for the job
#PBS -N ngmix_run
# Join output and errors in one file
#PBS -j oe
# Queue
## #PBS -q batch
# Request number of cores
#PBS -l nodes=1:ppn=1
#PBS -l walltime=48:00:00

# Activate conda environment
module load intelpython/3
source activate $HOME/.conda/envs/shapepipe

# Run ShapePipe
cd $HOME/ShapePipe
$HOME/.conda/envs/shapepipe/bin/python /home/guinot/ShapePipe/shapepipe/modules/ngmix_runner_script.py $PBS_ARRAYID /home/guinot/ShapePipe/example/test_full_auto/ngmix_process_list/process_list.npy /s03data2/guinot/pipeline_output/log_exp_headers.sqlite /home/guinot/ShapePipe_dir2/ngmix_output

# Return exit code
exit 0
