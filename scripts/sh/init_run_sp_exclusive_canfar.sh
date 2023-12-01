#!/bin/bash

echo "start init run sp exclusive canfar"

. /opt/conda/etc/profile.d/conda.sh

conda activate shapepipe

basedir=$HOME/cosmostat/P3_v2/psfex
cd $basedir


ID=$1
n_SMP=$2
typ=$3
echo "ID=$ID n_SMP=$n_SMP type=$typ"

export SP_RUN=.
export SP_CONFIG=$HOME/shapepipe/example/cfis
shapepipe_run -c $SP_CONFIG/config_exp_Pi.ini -e $ID

cd $basedir

echo "end init run sp exclusive canfar"

