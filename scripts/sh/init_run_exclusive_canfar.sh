#!/bin/bash

echo "start init run exclusive canfar"

#echo init_canfar > ~/init_canfar.log
#date >> ~/init_canfar.log

. /opt/conda/etc/profile.d/conda.sh

conda activate shapepipe

basedir=$HOME/cosmostat/P3_v2/psfex
cd $basedir


ID=$1
n_SMP=$2
typ=$3
echo "ID=$ID n_SMP=$n_SMP type=$typ"

cd ${typ}_runs

mkdir $ID
cd $ID
mkdir output
cd output
cp $basedir/output/log_run_sp.txt .
ln -s $basedir/output/log_exp_headers.sqlite
for dir in $basedir/output/run_sp_*; do
	ln -s $dir
done
cd ..

#export SP_RUN=.
#export SP_CONFIG=$HOME/shapepipe/example/cfis
#shapepipe_run -c $SP_CONFIG/config_tile_Sx.ini -e $ID

job_sp_canfar.bash -p psfex -j 32 -e $ID -n $n_SMP

cd $basedir

echo "end init run tile canfar"

