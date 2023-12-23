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
VERBOSE=1


# TODO: psf

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

if [ "$kind" == "-1" ]; then
  echo "No image kind indicated, use option -k"
  exit 4
fi

# Functions

## Print string, executes command, and prints return value.
function command () {
   cmd=$1
   dry_run=$2

   RED='\033[0;31m'
   GREEN='\033[0;32m'
   NC='\033[0m' # No Color
   # Color escape characters show up in log files
   #RED=''
   #GREEN=''
   #NC=''

   if [ $VERBOSE == 1 ]; then
        echo "running '$cmd' (dry run=$dry_run)"
   fi
   if [ "$dry_run" == "0" ]; then
        $cmd
        res=$?

        if [ $VERBOSE == 1 ]; then
            if [ $res == 0 ]; then
              echo -e "${GREEN}success, return value = $res${NC}"
            else
              echo -e "${RED}error, return value = $res${NC}"
              if [ $STOP == 1 ]; then
                  echo "${RED}exiting  $(basename "$0")', error in command '$cmd'${NC}"
                  exit $res
              else
                  echo "${RED}continuing '$(basename "$0")', error in command '$cmd'${NC}"
              fi
            fi
        fi
   fi
}

echo "start init_run_exclusive_canfar"

if [ "$dry_run" == 1 ]; then
  echo "in dry run mode"
fi

. /opt/conda/etc/profile.d/conda.sh

conda activate shapepipe

basedir=$HOME/cosmostat/P3_v2/psfex
cd $basedir


# Update links to exposure run directories
if [ "$kind" == "tile" ]; then
  command "rm -rf tile_runs/$ID/output/run_exp_SxSePsf*" $dry_run
  command "link_to_exp_for_tile.py -t $ID -i tile_runs -I exp_runs" $dry_run
fi


cd ${kind}_runs

if [ ! -d "$ID" ]; then
  command "mkdir $ID" $dry_run
fi

cd $ID
pwd

if [ ! -d "output" ]; then
  command "mkdir output" $dry_run
fi

cd output

if [ ! -f log_exp_headers.sqlite ]; then
  command "ln -s $basedir/output/log_exp_headers.sqlite" $dry_run
fi

# Remove potentially obsolete link
#rm run_sp_exp_SpMh*
#rm  run_sp_MaMa_*

# Update links to global run directories (GiFeGie, Uz, Ma?, combined_flag?)
for dir in $basedir/output/run_sp_*; do
  command "ln -sf $dir" $dry_run
done
if [ ! -e run_sp_combined_flag ]; then
  command "ln -s $basedir/output/run_sp_combined_flag" $dry_run
fi

(( do_job= $job & 128 ))
#if [[ $do_job != 0 ]]; then
if [ 0 == 1 ]; then

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
      command "rm -rf run_sp_tile_ngmix_Ng${k}u" $dry_run
    done
  else
    if [ "$job" == "128" ]; then
      echo "ngmix found complete, all is well, exiting"
      exit 0
    fi
  fi

fi 

cd ..
command update_runs_log_file.py $dry_run

echo -n "pwd: "
pwd

echo -n "environment: "
echo $CONDA_PREFIX

command "job_sp_canfar.bash -p psfex -j $job -e $ID --n_smp $N_SMP" $dry_run

cd $basedir

echo "end init run tile canfar"

