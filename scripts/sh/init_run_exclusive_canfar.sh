#!/bin/bash

# init_run_exclusive_canfar.sh

# Command line arguments
## Default values
job=-1
ID=-1
N_SMP=1
dry_run=0
nsh_jobs=8
dir=`pwd`

# mh_local is 0 (1) if merge_header_runner is run on all exposures,
# which is standard so far (run on exposures of given tile only; new)
mh_local=0
VERBOSE=1


## Help string
usage="Usage: $(basename "$0") -j JOB -e ID -k KIND [OPTIONS]
\n\nOptions:\n
   -h\tthis message\n
   -j, --job JOB\tRUnning JOB, bit-coded\n
   -e, --exclusive ID
    \timage ID\n
   -p, --psf MODEL\n
    \tPSF model, one in ['psfex'|'mccd'], default='$psf'\n
   -m, --mh_local MH\n
   \tmerged header file local (MH=0) or global (MH=1); default is $MH\n
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
    -p|--psf)                                                                   
      psf="$2"                                                                  
      shift                                                                     
      ;; 
    -m|--mh_local)                                                                   
      mh_local="$2"                                                                  
      shift                                                                     
      ;; 
    -N|--N_SMP)                                                                 
      N_SMP="$2"                                                                
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

## Check options
if [ "$job" == "-1" ]; then
  echo "No job indicated, use option -j"
  exit 2
fi

if [ "$exclusive" == "-1" ]; then
  echo "No image ID indicated, use option -e"
  exit 3
fi

if [ "$psf" != "psfex" ] && [ "$psf" != "mccd" ]; then
  echo "PSF (option -p) needs to be 'psfex' or 'mccd'"
  exit 4
fi

if [ "$mh_local" != "0" ] && [ "$mh_local" != "1" ]; then
  echo "MH (option -m) needs to be 0 or 1"
  exit 5
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

# Set kind
job_to_test=16
kind="none"

# loop over possible job numbers
while  [ $job_to_test -le 1024 ]; do

  (( do_job = $job & $job_to_test ))
  if [[ $do_job != 0 ]]; then
    
      if [ $job_to_test == 32 ]; then
        if [ "$kind" == "tile" ]; then
          echo "Error: Invalid job $job. mixing tile and exp kinds"
          exit 6
        fi

        # job=32 -> set kind to exp
        kind="exp"
      else
        if [ "$kind" == "exp" ]; then
          echo "Error: Invalid job $job. mixing tile and exp kinds"
          exit 6
        fi

        # job != 32 -> set kind to tile
        kind="tile"
      fi

  fi

  # Multiply job number by two to get next biwise number
  job_to_test=$((job_to_test * 2))
done

if [ "$kind" == "none" ]; then
  echo "Error: invalid job $job"
  exit 5
fi


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

if [ "$mh_local" == "0" ]; then
  if [ ! -f log_exp_headers.sqlite ]; then
    # Global Mh and file does not exist -> symlink to
    # gllobal mh file
    command "ln -s $dir/output/log_exp_headers.sqlite" $dry_run
  fi
fi


# Update links to global run directories (GiFeGie, Uz)
for my_dir in $dir/output/run_sp_[GU]*; do
  command "ln -sf $my_dir" $dry_run
done
# Combined flags
command "ln -sf $dir/output/run_sp_Ma_tile" $dry_run
command "ln -sf $dir/output/run_sp_Ma_exp" $dry_run
# exp Sp
command "ln -sf $dir/output/run_sp_exp_SpMh" $dry_run


(( do_job = $job & 16 ))
if [ "$mh_local" == "1" ] && [ $do_job != 0 ]; then

  # Remove previous Sx runs
  command "rm -rf $dir/output/run_tile_Sx_*" $dry_run

  if [ "$ID" == "-1" ]; then
    echo "ID needs to be given (option -e) for mh_local and job&16"
    exit 6 
  fi

  if [ -L log_exp_headers.sqlite ]; then
    # Local Mh and symlink -> remove previous link to
    # (potentially incomplete) global mh file
    echo "Removing previous mh sym link"
    command "rm log_exp_headers.sqlite" $dry_run
  else
    echo "no mh link found"
  fi

  echo "Creating local mh file"

  # Remove previous (local) split_exp dir
  command "rm -rf run_sp_exp_Sp" $dry_run

  # Create new split exp run dir
  new_dir="run_sp_exp_Sp//split_exp_runner/output"
  command "mkdir -p $new_dir" $dry_run

  # Link to all header files of exposures used for the current tile
  IDs=`echo $ID | tr "." "-"`
  for exp_ID in `cat run_sp_GitFeGie_*/find_exposures_runner/output/exp_numbers-$IDs.txt` ; do
    x=`echo $exp_ID | tr -d p `
    command "ln -s $dir/output/run_sp_exp_SpMh/split_exp_runner/output/headers-$x.npy $new_dir/headers-$x.npy" $dry_run
  done

  # Run merge_headers_runner on local exposure selection
  cd ..
  command "update_runs_log_file.py" $dry_run
  export SP_RUN=`pwd`
  command "shapepipe_run -c cfis/config_exp_Mh.ini" $dry_run
  cd output

  # Remove previous Sextractor run
  command "rm -rf run_sp_tile_Sx_*" $dry_run
fi

# Update links to exposure run directories, which were created in job 32
(( do_job = $job & 64 ))
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
      command "rm -rf `ls -rt1d run_sp_tile_Sx_* | head -$n_remove`" $dry_run
    fi

    # Remove previous runs of this job
    rm -rf run_sp_tile_PsViSmVi*
  fi
fi

(( do_job = $job & 256 ))
if [[ $do_job != 0 ]]; then

  # Remove previous runs of this job
  rm -rf run_sp_Ms_20??_*

fi

(( do_job = $job & 512 ))
if [[ $do_job != 0 ]]; then

  # Remove previous runs of this job
  rm -rf run_sp_Mc_20??_*

fi


cd ..

# Update log file
command update_runs_log_file.py $dry_run

echo -n "pwd: "
pwd

echo -n "environment: "
echo $CONDA_PREFIX

# To avoid (new?) qt error with setools (-j 32)
export DISPLAY=:1.0

command "job_sp_canfar.bash -p psfex -j $job -e $ID --n_smp $N_SMP" $dry_run

cd $dir

echo "end init run tile canfar"
