#!/bin/bash

# init_run_exclusive_canfar.sh

# Command line arguments
## Default values
job=-1
ID=-1
n_SMP=1
kind=-1

## Help string
usage="Usage: $(basename "$0") -j JOB -e ID  -k KIND [OPTIONS]
\n\nOptions:\n
   -h\tthis message\n
   -j, --job JOB\tRUnning JOB, bit-coded\n
   -e, --exclusive ID
    \timage ID\n
   -p, --psf MODEL\n
    \tPSF model, one in ['psfex'|'mccd'], default='$psf'\n
   -k, --kind KIND\n
    \timage kind, allowed are 'tile' and 'exp'\n
   -n, --n_SMP N_SMOp\n
    \tnumber of jobs (SMP mode only), default from original config files\n
"

## Help if no arguments                                                         
if [ -z $1 ]; then                                                              
        echo -ne $usage                                                         
        exit 1                                                                  
fi 

## Parse command line                                                           
while [ $# -gt 0 ]; do                                                          
  case "$1" in                                                                  
    -h)                                                                         
      echo -ne $usage                                                           
      exit 0                                                                    
      ;;                                                                        
    -j|--job)                                                                   
      job="$2"                                                                  
      shift                                                                     
      ;; 
    -e|--exclusive)                                                             
      ID="$2"                                                            
      shift                                                                     
      ;;
    -n|--n_SMP)                                                                 
      n_SMP="$2"                                                                
      shift                                                                     
      ;;                                                                        
    -k|--kind)                                                                 
      kind="$2"                                                                
      shift                                                                     
      ;;                                                                        
  esac                                                                          
  shift                                                                         
done

# Check options
if [ "$job" == "-1" ]; then
  echo "No job indicated, use option -j"
  exit 2
fi

if [ "$exclusive" == "-1" ]; then
  echo "No image ID indicated, use option -e"
  exit 3
fi

if [ "kind" == "-1" ]; then
  echo "No image kind indicated, use option -k"
  exit 4
fi

echo "start init run exclusive canfar"

. /opt/conda/etc/profile.d/conda.sh

conda activate shapepipe

basedir=$HOME/cosmostat/P3_v2/psfex
cd $basedir


echo "ID=$ID n_SMP=$n_SMP kind=$kind"

if [ "$kind" == "tile" ]; then
  rm -rf tile_runs/$ID/output/run_exp_SxSePsf*
  link_to_exp_for_tile.py -t $ID -i tile_runs -I exp_runs
fi

cd ${kind}_runs

if [ ! -d "$ID" ]; then
  mkdir $ID
fi

cd $ID

if [ ! -d "output" ]; then
  mkdir output
fi

cd output

if [ ! -f log_exp_headers.sqlite ]; then
  ln -s $basedir/output/log_exp_headers.sqlite
fi

# Remove potentially obsolete link
#rm run_sp_exp_SpMh*
#rm  run_sp_MaMa_*

for dir in $basedir/output/run_sp_*; do
 ln -sf $dir
done
if [ ! -e run_sp_combined_flag ]; then
  ln -s $basedir/output/run_sp_combined_flag
fi

cd ..
update_runs_log_file.py

echo -n "pwd: "
pwd

#export SP_RUN=.
#export SP_CONFIG=$HOME/shapepipe/example/cfis
#shapepipe_run -c $SP_CONFIG/config_exp_Pi.ini -e $ID

echo -n "environment: "
echo $CONDA_PREFIX

cmd="job_sp_canfar.bash -p psfex -j $job -e $ID -n $n_SMP"
echo "Running commnd $cmd ..."

$cmd

cd $basedir

echo "end init run tile canfar"

