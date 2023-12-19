#!/bin/bash

# init_run_exclusive_canfar.sh

# Command line arguments
## Default values
job=-1
ID=-1
N_SMP=1
kind=-1
dry_run=0
nsh_jobs=8

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
   -N, --N_SMP N_SMOp\n
    \tnumber of jobs (SMP mode only), default from original config files\n
   -n, --dry_run\n
    \tdry run, no actuall processing\n
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
    -N|--N_SMP)                                                                 
      n_SMP="$2"                                                                
      shift                                                                     
      ;;                                                                        
    -k|--kind)                                                                 
      kind="$2"                                                                
      shift                                                                     
      ;;                                                                        
    -n|--dry_run)
      dry_run=1
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

if [ "$dry_run" == 1 ]; then
  echo "in dry run mode"
fi

. /opt/conda/etc/profile.d/conda.sh

conda activate shapepipe

basedir=$HOME/cosmostat/P3_v2/psfex
cd $basedir


if [ "$dry_run" == "0" ]; then

  if [ "$kind" == "tile" ]; then
    rm -rf tile_runs/$ID/output/run_exp_SxSePsf*
    link_to_exp_for_tile.py -t $ID -i tile_runs -I exp_runs
  fi

  cd ${kind}_runs

  if [ ! -d "$ID" ]; then
    mkdir $ID
  fi

  cd $ID
  pwd

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

  # Indentify and remove unfinished ngmix dirs
  min_n_out=2
  for k in $(seq 1 $nsh_jobs); do
    ngmix_run="run_sp_tile_ngmix_Ng${k}u/ngmix_runner"
    if [ -e "$ngmix_run" ]; then
      ngmix_out="$ngmix_run/output"
      n_out=`ls -rlt $ngmix_out | wc -l`
      if [ "$n_out" -lt "$min_n_out" ]; then
          min_n_out=$n_out
      fi
      #echo $k $n_out $min_n_out
    else
      echo "ngmix separated run #$k not found"
      min_n_out=0
    fi
  done
  if [ "$min_n_out" -lt "2" ]; then
    echo "At least one ngmix separated run no output files"
    for k in $(seq 1 $nsh_jobs); do
      cmd="rm -rf run_sp_tile_ngmix_Ng${k}u"
      echo $cmd
      `$cmd`
    done
  else
    if [ "$job" == "128" ]; then
      echo "ngmix found complete, all is well, exiting"
      exit 0
    fi
  fi

  cd ..
  update_runs_log_file.py

  echo -n "pwd: "
  pwd

fi

#export SP_RUN=.
#export SP_CONFIG=$HOME/shapepipe/example/cfis
#shapepipe_run -c $SP_CONFIG/config_exp_Pi.ini -e $ID

echo -n "environment: "
echo $CONDA_PREFIX

cmd="job_sp_canfar.bash -p psfex -j $job -e $ID --n_smp $N_SMP"
echo -n "Running commnd '$cmd' ..."

if [ "$dry_run" == 1 ]; then
  echo " dry run"
else
  echo
  $cmd
fi

cd $basedir

echo "end init run tile canfar"

