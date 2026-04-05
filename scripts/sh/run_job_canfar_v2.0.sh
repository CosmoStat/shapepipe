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
check=0
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
   --check\t\tcheck download completeness only (job 8), no processing\n
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
    --check)
      check=1
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



# Initialise exposure work directory: create dirs, exp_numbers file, config symlink.
# The exp_numbers-000-000.txt file is created only once (skipped if already exists).
# Args: $1 = exp_id, $2 = exp_work_dir
function init_exp_work_dir() {
  local exp_id=$1
  local exp_work_dir=$2
  local fe_output="$exp_work_dir/output/run_sp_Fe/find_exposures_runner/output"

  [ ! -d "$fe_output" ] && command "mkdir -p $fe_output" $dry_run
  [ ! -d "$exp_work_dir/output" ] && command "mkdir -p $exp_work_dir/output" $dry_run

  local exp_numbers_out="$fe_output/exp_numbers-000-000.txt"
  if [ ! -e "$exp_numbers_out" ] && [ "$dry_run" == "0" ]; then
    echo "$exp_id" > "$exp_numbers_out"
  fi

  if [ ! -e "$exp_work_dir/cfis" ]; then
    ln -sf ~/shapepipe/example/cfis "$exp_work_dir/cfis"
  fi
}


# Run a per-exposure job (e.g. job 8, 16).
# Args: $1 = job number, $2 = run_sp output dir prefix (e.g. "Gie")
#       $3 = space-separated list of "runner_subdir:N" completeness checks
#            all pairs must pass for an exposure to be considered complete
function run_job_exp() {
  local exp_job=$1
  local run_prefix=$2
  local complete_checks=$3

  exp_numbers_file=$(ls -t "$work_dir/output/run_sp_Fe"*/find_exposures_runner/output/"exp_numbers-${IDra}-${IDdec}.txt" 2>/dev/null | head -1)

  if [ -z "$exp_numbers_file" ]; then
    message "Exposure numbers file exp_numbers-${IDra}-${IDdec}.txt not found in $work_dir/output" $debug_out 10
  fi

  if [ "$check" == "1" ]; then
    message "Check mode: skipping job $exp_job" $debug_out -1
  fi

  local n_total=0
  local n_complete=0
  local n_incomplete=0

  # Loop over exposure IDs
  while IFS= read -r exp_id || [ -n "$exp_id" ]; do
    [ -z "$exp_id" ] && continue

    (( n_total++ ))

    # exp_id e.g. "2182795p": ab = first 2 chars, abcdefg = all but last char
    local exp_prefix="${exp_id:0:2}"
    local exp_base="${exp_id%?}"
    local exp_work_dir="$HOME/v${version}/exp/$exp_prefix/$exp_base"
    local exp_log_file="$exp_work_dir/job_sp_canfar_v2.0.log"

    # Create exp_numbers-000-000.txt and cfis link if not existent
    init_exp_work_dir "$exp_id" "$exp_work_dir"

    # Check completeness of existing run output
    local run_dir=$(ls -dt "$exp_work_dir/output/run_sp_${run_prefix}"* 2>/dev/null | head -1)
    local is_complete=1
    local check_desc=""
    for check_pair in $complete_checks; do
      local subdir="${check_pair%:*}"
      local n_threshold="${check_pair##*:}"
      local n_out=0

      out_dir="$run_dir/${subdir}/output"

      # Remove broken symlinks in module output dir
      for f in "$out_dir"/*; do
        if [ -L "$f" ] && [ ! -e "$f" ]; then
          message "Removing broken link: $f" $debug_out -1
          command "rm $f" $dry_run
        fi
      done

      [ -n "$run_dir" ] && n_out=$(ls "$out_dir/" 2>/dev/null | wc -l)
      check_desc+="${subdir}:${n_out}/${n_threshold} "
      [ "$n_out" -lt "$n_threshold" ] && is_complete=0
    done

    if [ "$is_complete" == "1" ]; then
      message "Skipping $exp_id: run_sp_${run_prefix} complete ( $check_desc)" $debug_out -1
      (( n_complete++ ))
      continue
    fi

    # Report incomplete/missing in check mode; in run mode handle and proceed
    if [ "$check" == "1" ]; then
      if [ -n "$run_dir" ]; then
        message "  incomplete: $exp_id ($check_desc)" $debug_out -1
      else
        message "  missing: $exp_id" $debug_out -1
      fi
      (( n_incomplete++ ))
      continue
    fi

    cd "$exp_work_dir"
    command "update_runs_log_file.py" $dry_run

    echo "$(basename "$0") -j $exp_job -e $exp_id" > "$exp_log_file"
    command "job_sp_canfar_v2.0.bash -p $psf -j $exp_job --n_smp $N_SMP --nsh_jobs $N_SMP --debug_out $debug_out" $dry_run 2>&1 | tee -a "$exp_log_file"

    cd "$dir"

  done < "$exp_numbers_file"

  message "Tile $ID job $exp_job: $n_complete/$n_total exposures complete" $debug_out -1
  if [ "$n_incomplete" -gt 0 ]; then
    message "WARNING: $n_incomplete/$n_total exposures incomplete or missing" $debug_out -1
  fi
}


# Init message
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
[ ! -d "$work_dir" ] && command "mkdir -p $work_dir" $dry_run
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

echo -n "pwd: "; pwd

# Avoid Qt error with setools
export DISPLAY=:1.0


# Run job — on scratch if available, otherwise in-place
# Job 1: download tile images (config_Git_vos.ini)

if [ "$scratch" != "-1" ]; then
  scratch_work="$scratch/tiles/$IDra/$ID"
  message "Copying work dir to scratch: $scratch_work" $debug_out -1
  [ ! -d "$scratch/tiles/$IDra" ] && command "mkdir -p $scratch/tiles/$IDra" $dry_run
  command "cp -aR . $scratch_work" $dry_run
  command "cd $scratch_work" $dry_run
fi

IDdec=${ID##*.}

if [ "$job" == "8" ]; then

  # Job 8: retrieve exposure images (config_Gie_vos.ini)
  run_job_exp $job "Gie" "get_images_runner:6"

elif [ "$job" == "16" ]; then

  # Job 16: split images and merge headers (config_exp_SpMh.ini)
  run_job_exp $job "exp_SpMh" "split_exp_runner:120 merge_headers_runner:1"

elif [ "$job" == "32" ]; then

  # Job 32: mask exposures (config_exp_Ma_onthefly.ini)
  run_job_exp $job "exp_Ma" "mask_runner:4"

else

  echo "$(basename "$0") $@" > "$log_file"
  command "job_sp_canfar_v2.0.bash -p $psf -j $job --n_smp $N_SMP --nsh_jobs $N_SMP --debug_out $debug_out" $dry_run 2>&1 | tee -a "$log_file"

fi

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
