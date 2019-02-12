#!/bin/bash
#$ -M <name>@cea.fr
#$ -m ea
#$ -N shapepipe_run
#$ -P P_euclid
#$ -j y
#$ -l os=cl7
#$ -pe multicores 4

ccenv anaconda
source activate $HOME/.conda/envs/shapepipe

cd $HOME/ShapePipe

mpiexec -n $NSLOTS $HOME/.conda/envs/shapepipe/bin/python shapepipe_run.py -c example/config_mpi.ini

exit 0
