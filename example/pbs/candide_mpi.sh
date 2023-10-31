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
export SPENV="$HOME/.conda/envs/shapepipe"

# Full path to example config file and input data
export SPDIR="$HOME/shapepipe"

# Load modules
module remove gcc
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
  NSLOTS=8
  echo "Using $NSLOTS CPUs set by hand"
fi

# Test: print hostname.
# Only version 5.0.0 downloaded from the web recognised the --mca argument
#/home/mkilbing/bin/mpirun -np $NSLOTS --mca pml ob1 --mca btl ^openib --mca orte_base_help_aggregate 0 --mca plm_tm_verbose 1 hostname

#/softs/openmpi/5.0.0-torque-CentOS7/bin/mpirun -map-by node --mca pml ob1 --mca btl ^openib --mca orte_base_help_aggregate 0 --mca plm_tm_verbose 1 hostname
#/softs/openmpi/5.0.0-torque-CentOS7/bin/mpirun -map-by node $SPENV/bin/shapepipe_run -c $SPDIR/example/pbs/config_mpi.ini
#/home/mkilbing/bin/mpirun -map-by node $SPENV/bin/shapepipe_run -c $SPDIR/example/pbs/config_mpi.ini
#/home/mkilbing/bin/mpirun -n $NSLOTS --mca pml ob1 --mca btl ^openib --mca orte_base_help_aggregate 0 --mca plm_tm_verbose 1 $SPENV/bin/shapepipe_run -c $SPDIR/example/pbs/config_mpi.ini
#$SPENV/bin/mpiexec -n $NSLOTS $SPENV/bin/shapepipe_run -c $SPDIR/example/pbs/config_mpi.ini
#/softs/openmpi/5.0.0-torque-CentOS7/bin/mpirun -np $NSLOTS --mca pml ob1 --mca btl ^openib --mca orte_base_help_aggregate 0 --mca plm_tm_verbose 1 $SPENV/bin/shapepipe_run -c $SPDIR/example/pbs/config_mpi.ini
/softs/openmpi/5.0.0-torque-CentOS7/bin/mpirun -np $NSLOTS hostname
/softs/openmpi/5.0.0-torque-CentOS7/bin/mpirun -np $NSLOTS $SPENV/bin/shapepipe_run -c $SPDIR/example/pbs/config_mpi.ini

#/home/mkilbing/bin/mpirun -np $NSLOTS --mca pml ob1 --mca btl ^openib --mca orte_base_help_aggregate 0 --mca plm_tm_verbose 1 $SPENV/bin/shapepipe_run -c $SPDIR/example/pbs/config_mpi.ini

# Return exit code
exit 0
