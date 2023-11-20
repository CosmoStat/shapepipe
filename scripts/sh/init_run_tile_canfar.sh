#!/bin/bash

echo "start init run tile canfar"

#echo init_canfar > ~/init_canfar.log
#date >> ~/init_canfar.log

. /opt/conda/etc/profile.d/conda.sh

conda activate shapepipe

basedir=cosmostat/P3_v2/psfex 
cd $basedir

tile_ID=$1
n_SMP=$2
echo "tile_ID=$tile_ID n_SMP=$n_SMP"

mkdir $tile_ID
cd $tile_ID
mkdir output
cd output
cp $basedir/log_run_sp.txt .
for dir in $basedir/run_sp_*; do
	ln -s $$dir
done
cd ..

job_sp $tile_ID -p psfex -j 16 -n $n_SMP

echo "end init run tile canfar"

