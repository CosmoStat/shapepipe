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
dir=`pwd`
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
   -d, --directory\n
    \trun directory, default is pwd ($dir)\n
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
      N_SMP="$2"                                                                
      shift                                                                     
      ;;                                                                        
    -k|--kind)                                                                 
      kind="$2"                                                                
      shift                                                                     
      ;;                                                                        
    -d|--directory)
      dir="$2"
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

cd $dir
echo $pwd

if [ ! -d ${kind}_runs ]; then
  command "mkdir ${kind}_runs" $dry_run
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
  command "ln -s $dir/output/log_exp_headers.sqlite" $dry_run
fi


# Update links to global run directories (GiFeGie, Uz, Ma?, combined_flag?)
for dir in $dir/output/run_sp_*; do
  command "ln -sf $dir" $dry_run
done

# Update links to exposure run directories, which were created in job 32
(( do_job= $job & 64 ))
if [[ $do_job != 0 ]]; then
  if [ "$kind" == "tile" ]; then
    cd ../../..
    command "link_to_exp_for_tile.py -t $ID -i tile_runs -I exp_runs" $dry_run
    cd ${kind}_runs/$ID/output

    # Remove duplicate job-32 runs (tile detection)
    n_32=`ls -rt1d run_sp_tile_Sx_* | wc -l`
    if [ "$n_32" != "1" ]; then
      n_remove="$(($n_32-1))"
      echo "removing $n_remove duplicate old job-32 runs"
      rm -rf `ls -rt1d run_sp_tile_Sx_* | head -$n_remove`
    fi

    # Remove previous runs of this job
    rm -rf run_sp_tile_PsViSmVi*
  fi
fi

(( do_job= $job & 256 ))
if [[ $do_job != 0 ]]; then

  # Remove previous runs of this job
  rm -rf run_sp_Ms_20??_*
  rm -rf run_sp_Mc_20??_*

fi

cd ..

# Update log file
command update_runs_log_file.py $dry_run

echo -n "pwd: "
pwd

echo -n "environment: "
echo $CONDA_PREFIX

command "job_sp_canfar.bash -p psfex -j $job -e $ID --n_smp $N_SMP" $dry_run

cd $dir

echo "end init run tile canfar"
