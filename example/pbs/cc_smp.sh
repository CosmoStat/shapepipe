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

$HOME/.conda/envs/sshapepipe/bin/python shapepipe_run.py -c example/pbs/config_smp.ini

exit 0
