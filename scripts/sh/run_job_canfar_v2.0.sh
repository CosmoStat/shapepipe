#!/bin/bash

# run_job_canfar_v2.0.sh
# Description: Initialise tile/exposure run directory and launch ShapePipe job.

# ShapePipe version
version="2.0"

# Command line arguments
## Default values
job=-1
ID=-1
psf='psfex'
N_SMP=1
dry_run=0
dir=`pwd`
debug_out=-1
#scratch="/scratch/$USER/shapepipe/v${version}"
scratch="-1"
test_only=0
VERBOSE=1

pat="-- "


## Help string
usage="Usage: $(basename "$0") -j JOB -e ID [OPTIONS]
\n\nOptions:\n
   -h\t\t\tthis message\n
   -j, --job JOB\t\trunning JOB, bit-coded\n
   -e, --exclusive ID\timage ID\n
   -p, --psf MODEL\n
   \t\t\tPSF model, one in ['psfex'|'mccd'], default='$psf'\n
   -N, --N_SMP N_SMP\tnumber of SMP jobs, default from original config files\n
   -d, --directory DIR\trun directory, default is pwd ($dir)\n
   -S, --scratch DIR\tprocessing scratch directory, default='$scratch'; use -1 to disable\n
   -n, --dry_run\t\tdry run, no actual processing; default is $dry_run\n
   --debug_out PATH\tdebug output file PATH, default not used\n
   --test\t\ttest mode, no processing\n
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
    -N|--N_SMP)
      N_SMP="$2"
      shift
      ;;
    -d|--directory)
      dir="$2"
      shift
      ;;
    -S|--scratch)
      scratch="$2"
      shift
      ;;
    -n|--dry_run)
      dry_run="$2"
      shift
      ;;
    --debug_out)
      debug_out="$2"
      shift
      ;;
    --test)
      test_only=1
      ;;
  esac
  shift
done


# Functions

function message() {
  msg=$1
  my_debug_out=$2
  my_exit=$3

  echo $msg
  if [ "$my_debug_out" != "-1" ]; then
    echo ${pat}$msg >> $my_debug_out
  fi

  if [ "$my_exit" != "-1" ]; then
    echo "${pat}exiting with code $my_exit" >> $my_debug_out
    exit $my_exit
  fi
}


# Init message
message "test=$test_only" $debug_out -1
if [ "$test_only" == "1" ]; then
  message "$(basename "$0") test mode, exiting." $debug_out 0
else
  message "$(basename "$0") processing mode, starting." $debug_out -1
fi


## Check options
message "checking options" $debug_out -1

if [ "$job" == "-1" ]; then
  message "No job indicated, use option -j" $debug_out 2
fi

if [ "$ID" == "-1" ]; then
  message "No image ID indicated, use option -e" $debug_out 3
fi

if [ "$psf" != "psfex" ] && [ "$psf" != "mccd" ]; then
  message "PSF (option -p) needs to be 'psfex' or 'mccd'" $debug_out 4
fi

if [ "$dry_run" != "0" ] && [ "$dry_run" != "1" ]; then
  message "dry_run must be 0 or 1, not $dry_run" $debug_out 8
fi


# Start script

source $HOME/shapepipe/scripts/sh/functions.sh

message "Starting $(basename "$0")" $debug_out -1
message "`date`" $debug_out -1
message "ID=$ID" $debug_out -1

if [ "$dry_run" == "1" ]; then
  message "running in dry run mode" $debug_out -1
fi

CONDA_PREFIX=$HOME/.conda/envs/shapepipe
PATH=$PATH:$CONDA_PREFIX/bin

cd $dir

# Derive tile path components from ID (e.g. "000.227" -> IDra="000")
IDra=${ID%%.*}
work_dir="$dir/tiles/$IDra/$ID"
log_file="$work_dir/job_sp_canfar_v2.0.log"

# Create tile work directory
command "mkdir -p $work_dir" $dry_run
cd $work_dir

# Config symlink
ln -sf ~/shapepipe/example/cfis

# Write ID to first input
if [ ! -e tile_numbers.txt ]; then
  echo $ID > tile_numbers.txt
fi

# Output directory
if [ ! -d "output" ]; then
  command "mkdir output" $dry_run
fi

# Update log file
#command "update_runs_log_file.py" $dry_run

echo -n "pwd: "; pwd

# Avoid Qt error with setools
export DISPLAY=:1.0


# Run job — on scratch if available, otherwise in-place
# Job 1: download tile images (config_Git_vos.ini)

if [ "$scratch" != "-1" ]; then
  scratch_work="$scratch/tiles/$IDra/$ID"
  message "Copying work dir to scratch: $scratch_work" $debug_out -1
  command "mkdir -p $scratch/tiles/$IDra" $dry_run
  command "cp -aR . $scratch_work" $dry_run
  command "cd $scratch_work" $dry_run
fi

echo "$(basename "$0") $@" > "$log_file"
command "job_sp_canfar_v2.0.bash -p $psf -j $job --n_smp $N_SMP --nsh_jobs $N_SMP --debug_out $debug_out" $dry_run 2>&1 | tee -a "$log_file"

if [ "$scratch" != "-1" ]; then
  message "Syncing output from scratch back to permanent dir" $debug_out -1
  command "rsync -a output/ $work_dir/output/" $dry_run
  command "cd $work_dir" $dry_run
  command "update_runs_log_file.py" $dry_run
  command "rm -rf $scratch_work" $dry_run
  command "cd $dir" $dry_run
fi

cd $dir

message "End $(basename "$0")" $debug_out -1
