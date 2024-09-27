#!/bin/bash

##########################
# MPI Script for CANDIDE #
##########################

# Receive email when job finishes or aborts
## #PBS -M <name>@cea.fr
## #PBS -m ea

# Set a name for the job
#PBS -N shapepipe_mpi

# Join output and errors in one file
#PBS -j oe

# Set maximum computing time (e.g. 5min)
#PBS -l walltime=00:05:00

# Request number of cores (e.g. 2 from 2 different machines)
#PBS -l nodes=2:ppn=2

# Full path to environment
export SPENV="$HOME/.conda/envs/shapepipe_mpi"

# Full path to example config file and input data
export SPDIR="$HOME/shapepipe"

# Load modules
module load gcc/9.3.0
module load intelpython/3-2023.1.0
module load openmpi/5.0.0

# Activate conda environment
source activate $SPENV

# Other options to test
# -map-by

if [ -f "$PBS_NODEFILE" ]; then
  NSLOTS=`cat $PBS_NODEFILE | wc -l`
  echo "Using $NSLOTS CPUs from PBS_NODEFILE $PBS_NODEFILE"
else
  NSLOTS=4
  echo "Using $NSLOTS CPUs set by hand"
fi

# Creates #node output dirs
MPI_CMD=/softs/openmpi/5.0.0-torque-CentOS7/bin/mpirun
MPI_ARGS="-np $NSLOTS"

${MPI_CMD} ${MPI_ARGS} hostname
${MPI_CMD} ${MPI_ARGS} $SPENV/bin/shapepipe_run -c $SPDIR/example/pbs/config_mpi.ini

# Return exit code
exit 0
