#!/bin/bash

echo "start init canfar"

#echo init_canfar > ~/init_canfar.log
#date >> ~/init_canfar.log

. /opt/conda/etc/profile.d/conda.sh

conda activate shapepipe

cd cosmostat/P3_v2/psfex 

tile_ID=$1
n_SMP=$2
echo "tile_ID=$tile_ID n_SMP=$n_SMP"

job_sp $tile_ID -p psfex -j 8 -n $n_SMP

echo "end init canfar"

