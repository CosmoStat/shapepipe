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

#mkdir $ID

cd $ID

if [ ! -d "output" ]; then
  mkdir output
fi

  cd output
  #ln -s $basedir/output/log_exp_headers.sqlite
  # Remove potentially obsolete link
  rm run_sp_exp_SpMh*
  for dir in $basedir/output/run_sp_*; do
	  ln -s $dir
  done
  rm  run_sp_MaMa_*
  ln -s $basedir/output/run_sp_combined_flag
  cd ..
  update_runs_log_file.py

  pwd

#export SP_RUN=.
#export SP_CONFIG=$HOME/shapepipe/example/cfis
#shapepipe_run -c $SP_CONFIG/config_exp_Pi.ini -e $ID

job_sp_canfar.bash -p psfex -j 32 -e $ID -n $n_SMP

cd $basedir

echo "end init run tile canfar"

