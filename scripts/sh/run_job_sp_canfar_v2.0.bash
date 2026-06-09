#!/bin/bash

# run_job_sp_canfar_v2.0
# Description: Initialise tile/exposure run directory and launch ShapePipe job.

# Shared job-list description
source $HOME/shapepipe/scripts/sh/job_list_help.bash

# ShapePipe version
version="2.0"

# Command line arguments
## Default values
job=-1
ID=-1
psf='psfex'
tile_det='uc'
tile_mask=0
N_SMP=-1
dry_run=0
dir=`pwd`
debug_out=""

# Input type: data or image_sims
type="data"

#scratch="/scratch/$USER/shapepipe/v${version}"
scratch=""
test_only=0
check=0
force=0
VERBOSE=1

pat="-- "


## Help string
usage="Usage: $(basename "$0") -j JOB -e ID [OPTIONS]
\n\nOptions:\n
   -h\t\t\tthis message\n
   -j, --job JOB\t\trunning JOB, bit-coded\n
${JOB_LIST_HELP}   -e, --exclusive ID\timage ID\n
   -p, --psf MODEL\n
   \t\t\tPSF model, one in ['psfex'|'mccd'], default='$psf'\n
   --tile_det DET\t\ttile detection mode, one in ['sx'|'uc'], default='$tile_det'\n
   --tile_mask MASK\ttile masking, default='$tile_mask'\n
   -t, --type TYPE input type, allowed are 'data', 'image_sims', default='$type'\n
   -N, --N_SMP N_SMP\tnumber of SMP jobs, default from original config files\n
   -d, --directory DIR\trun directory, default is pwd ($dir)\n
   -S, --scratch DIR\tprocessing scratch directory, default=none\n
   -n, --dry_run\t\tDRY RUN, no actual processing; default is $dry_run\n
   --debug_out PATH\tdebug output file PATH, default=none\n
   --test\t\t\ttest mode, no processing\n
   --check\t\tcheck download completeness only (job 8), no processing\n
   --force\t\tremove existing module output dir(s) before running\n
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
    -t|--type)
      type="$2"
      shift
      ;;
    --tile_det)
      tile_det="$2"
      shift
      ;;
    --tile_mask)
      tile_mask="$2"
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
    --force)
      force=1
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
  if [ -n "$my_debug_out" ]; then
    echo ${pat}$msg >> $my_debug_out
  fi
  if [ -n "$log_file" ]; then
    echo ${pat}$msg >> $log_file
  fi

  if [ "$my_exit" != "-1" ]; then
    if [ -n "$my_debug_out" ]; then
      echo "${pat}exiting with code $my_exit" >> $my_debug_out
    else
      echo "${pat}exiting with code $my_exit"
    fi
    if [ -n "$log_file" ]; then
      echo "${pat}exiting with code $my_exit" >> $log_file
    fi
    exit $my_exit
  fi
}


# Initialise exposure work directory: create dirs, exp_numbers file, config symlink.
# The exp_numbers-000-000.txt file is created only once (skipped if already exists).
# Args: $1 = exp_id, $2 = exp_work_dir
function init_exp_work_dir() {
  local exp_id=$1
  local exp_work_dir=$2
  local fe_output="$exp_work_dir/output/run_sp_tile_Fe/find_exposures_runner/output"

  [ ! -d "$fe_output" ] && command "mkdir -p $fe_output" $dry_run
  [ ! -d "$exp_work_dir/output" ] && command "mkdir -p $exp_work_dir/output" $dry_run

  local exp_numbers_out="$fe_output/exp_numbers-000-000.txt"
  if [ ! -e "$exp_numbers_out" ] && [ "$dry_run" == "0" ]; then
    echo "$exp_id" > "$exp_numbers_out"
  fi

  if [ ! -e "$exp_work_dir/cfis" ]; then
    ln -sf $config_dir "$exp_work_dir/cfis"
  fi
}


# Run a per-exposure job (e.g. job 8, 16).
# Args: $1 = job number
#       $2 = space-separated list of run_sp_exp output dir prefixes (e.g. "Gie")
#            first prefix is the main one; all are force-removed when --force is set
#       $3 = space-separated list of completeness checks, each in one of two forms:
#              "runner_subdir:N[:subpath[:warn]]"
#                 check runner_subdir/output in the main (first) prefix run dir
#              "run_prefix:runner_subdir:N[:subpath[:warn]]"
#                 check runner_subdir/output in run_sp_exp_run_prefix* dir
#            all non-warn checks must pass for an exposure to be considered complete
function run_exp_job() {
  local exp_job=$1
  local run_prefixes=$2
  local complete_checks=$3
  local main_prefix="${run_prefixes%% *}"

  exp_numbers_file=$(ls -t "$work_dir/output/run_sp_tile_Fe"*/find_exposures_runner/output/"exp_numbers-${IDra}-${IDdec}.txt" 2>/dev/null | head -1)

  if [ -z "$exp_numbers_file" ]; then
    message "Exposure numbers file exp_numbers-${IDra}-${IDdec}.txt not found in $work_dir/output" "$debug_out" 10
  fi

  if [ "$check" == "1" ]; then
    message "Check mode: skipping job $exp_job" "$debug_out" -1
  fi

  local n_total=0
  local n_complete=0
  local n_incomplete=0
  local t_loop_start=$(date +%s)

  # Loop over exposure IDs
  while IFS= read -r exp_id || [ -n "$exp_id" ]; do
    [ -z "$exp_id" ] && continue

    (( n_total++ ))

    # exp_id e.g. "2182795p" (data) or "208659" (image_sims)
    # Strip trailing letter if present (data format); keep full id if numeric only.
    local exp_prefix="${exp_id:0:2}"
    local exp_base
    if [[ "${exp_id: -1}" =~ [a-zA-Z] ]]; then
      exp_base="${exp_id%?}"
    else
      exp_base="$exp_id"
    fi
    local exp_id_disp="${exp_prefix}/${exp_base}"
    local exp_work_dir="$dir/exp/$exp_prefix/$exp_base"
    local exp_log_file="$exp_work_dir/job_sp_canfar_v2.0.log"

    # Create exp_numbers-000-000.txt and cfis link if not existent
    init_exp_work_dir "$exp_id" "$exp_work_dir"

    # force: remove all existing run directories for each prefix before running
    if [ "$force" == "1" ]; then
      local run_prefix
      for run_prefix in $run_prefixes; do
        local dirs_to_remove
        dirs_to_remove=$(ls -d "$exp_work_dir/output/run_sp_exp_${run_prefix}"* 2>/dev/null)
        if [ -n "$dirs_to_remove" ]; then
          for d in $dirs_to_remove; do
            message "Force-removing $d" "$debug_out" -1
            command "rm -rf $d" $dry_run
          done
        fi
      done
    fi

    # Check completeness of existing run output (main prefix)
    local run_dir=$(ls -dt "$exp_work_dir/output/run_sp_exp_${main_prefix}"* 2>/dev/null | head -1)
    local is_complete=1
    local check_desc=""
    for check_pair in $complete_checks; do
      local field1="${check_pair%%:*}"
      local rest1="${check_pair#*:}"
      local field2="${rest1%%:*}"

      local check_run_dir subdir rest
      if [[ "$field2" =~ ^[0-9]+$ ]]; then
        # "subdir:N[...]" — check in main run dir
        check_run_dir="$run_dir"
        subdir="$field1"
        rest="$rest1"
      else
        # "run_prefix:subdir:N[...]" — check in that prefix's run dir
        check_run_dir=$(ls -dt "$exp_work_dir/output/run_sp_exp_${field1}"* 2>/dev/null | head -1)
        subdir="$field2"
        rest="${rest1#*:}"
      fi

      local n_threshold="${rest%%:*}"
      local rest2="${rest#*:}"
      local subpath=""
      local warn_only=0
      if [[ "$rest" == *:* ]]; then
        subpath="${rest2%%:*}"
        if [[ "$rest2" == *:* ]]; then
          [ "${rest2#*:}" == "warn" ] && warn_only=1
        fi
      fi
      local n_out=0

      local out_dir
      if [ -n "$subpath" ]; then
        out_dir="${check_run_dir}/${subdir}/output/${subpath}"
      else
        out_dir="${check_run_dir}/${subdir}/output"
      fi

      # Remove broken symlinks in module output dir
      for f in "$out_dir"/*; do
        if [ -L "$f" ] && [ ! -e "$f" ]; then
          message "Removing broken link: $f" "$debug_out" -1
          command "rm $f" $dry_run
        fi
      done

      [ -n "$check_run_dir" ] && n_out=$(ls "$out_dir/" 2>/dev/null | wc -l)
      check_desc+="${subdir}:${n_out}/${n_threshold} "
      if [ "$n_out" -lt "$n_threshold" ]; then
        if [ "$warn_only" == "1" ]; then
          #message "WARNING: ${subdir}: only $n_out/$n_threshold files" "$debug_out" -1
          : # do nothing
        else
          is_complete=0
        fi
      fi
    done

    if [ "$is_complete" == "1" ]; then
      message "Complete $exp_id_disp: run_sp_exp_${main_prefix} ( $check_desc)" "$debug_out" -1
      (( n_complete++ ))
      continue
    fi

    # Report incomplete/missing in check mode; in run mode handle and proceed
    if [ "$check" == "1" ]; then
      if [ -n "$run_dir" ]; then
        message "  Benign incomplete: $exp_id_disp ($check_desc)" "$debug_out" -1
      else
        message "  missing: $exp_id_disp" "$debug_out" -1
      fi
      (( n_incomplete++ ))
      continue
    fi

    cd "$exp_work_dir"
    command "update_runs_log_file.py" $dry_run

    # Build optional --debug_out flag (omit entirely when debug_out is empty)
    local debug_flag=""
    [ -n "$debug_out" ] && debug_flag="--debug_out $debug_out"

    echo "$(basename "$0") -j $exp_job -e $exp_id" > "$exp_log_file"
    echo "pwd=`pwd`"
    command "job_sp_canfar_v2.0.bash -c $config_dir -p $psf -r $retrieve --tile_det $tile_det --tile_mask $tile_mask -j $exp_job --n_smp $N_SMP --nsh_jobs $N_SMP $debug_flag" $dry_run 2>&1 | tee -a "$exp_log_file"
    echo "Done with job_sp_canfar_v2.0.bash"

  done < "$exp_numbers_file"

  local t_loop_end=$(date +%s)
  message "Exposure loop: $((t_loop_end - t_loop_start))s for $n_total exposures (job $exp_job)" "$debug_out" -1
  message "Tile $ID job $exp_job: $n_complete/$n_total exposures complete" "$debug_out" -1
  if [ "$n_incomplete" -gt 0 ]; then
    message "WARNING: $n_incomplete/$n_total exposures incomplete or missing" "$debug_out" -1
  fi
}


# Run a tile-level job (jobs 1, 2, 4, 128, 256, ...).
# Args: $1 = job number
#       $2 = space-separated list of run_sp_tile output dir prefixes (e.g. "Fe")
#            first prefix is the main one; all are force-removed when --force is set
#       $3 = space-separated completeness checks, each in one of two forms:
#              "runner_subdir:N[:subpath[:warn]]"
#                 check runner_subdir/output in the main (first) prefix run dir
#              "run_prefix:runner_subdir:N[:subpath[:warn]]"
#                 check runner_subdir/output in run_sp_tile_run_prefix* dir
#            omit or pass "" to skip the completeness check and always run

function run_tile_job() {
  local tile_job=$1
  local run_prefixes=$2
  local complete_checks=$3
  local main_prefix="${run_prefixes%% *}"

  # force: remove all existing run directories for each prefix before running
  if [ "$force" == "1" ]; then
    local run_prefix
    for run_prefix in $run_prefixes; do
      local dirs_to_remove
      dirs_to_remove=$(ls -d "$work_dir/output/run_sp_tile_${run_prefix}"* 2>/dev/null)
      if [ -n "$dirs_to_remove" ]; then
        for d in $dirs_to_remove; do
          message "Force-removing $d" "$debug_out" -1
          command "rm -rf $d" $dry_run
        done
      fi
    done
  fi

  # Locate most recent existing run directory for the main prefix
  local run_dir
  run_dir=$(ls -dt "$work_dir/output/run_sp_tile_${main_prefix}"* 2>/dev/null | head -1)

  local is_complete=1
  local check_desc=""

  if [ -n "$complete_checks" ]; then

    # turn inputs into arrays
    read -r -a run_prefix_array <<< "$run_prefixes"
    read -r -a check_array <<< "$complete_checks"

    # sanity check (optional but strongly recommended)
    if [ "${#run_prefix_array[@]}" -ne "${#check_array[@]}" ]; then
      echo "Error: run_prefixes and complete_checks must have same length"
      return 1
    fi

    # zip loop
    for i in "${!run_prefix_array[@]}"; do
      run_prefix="${run_prefix_array[$i]}"
      check_pair="${check_array[$i]}"

      local field1="${check_pair%%:*}"
      local rest1="${check_pair#*:}"
      local field2="${rest1%%:*}"

      local check_run_dir subdir rest
      check_run_dir=$(ls -dt "$work_dir/output/run_sp_tile_${run_prefix}"* 2>/dev/null | head -1)
      if [[ "$field2" =~ ^[0-9]+$ ]]; then
        subdir="$field1"
        rest="$rest1"
      else
        subdir="$field2"
        rest="${rest1#*:}"
      fi

      local n_threshold="${rest%%:*}"
      local rest2="${rest#*:}"
      local subpath=""
      local warn_only=0
      if [[ "$rest" == *:* ]]; then
        subpath="${rest2%%:*}"
        [[ "$rest2" == *:* ]] && [ "${rest2#*:}" == "warn" ] && warn_only=1
      fi

      local out_dir
      if [ -n "$subpath" ]; then
        out_dir="${check_run_dir}/${subdir}/output/${subpath}"
      else
        out_dir="${check_run_dir}/${subdir}/output"
      fi

      local n_out=0
      [ -n "$check_run_dir" ] && n_out=$(ls "$out_dir/" 2>/dev/null | wc -l)
      check_desc+="${run_prefix}/${subdir}:${n_out}/${n_threshold} "
      if [ "$n_out" -lt "$n_threshold" ]; then
        [ "$warn_only" != "1" ] && is_complete=0
      fi
    done
  fi

  if [ "$is_complete" == "1" ] && [ -n "$complete_checks" ]; then
    message "Complete: ( $check_desc)" "$debug_out" -1
    return 0
  fi

  if [ "$check" == "1" ]; then
    if [ -n "$run_dir" ]; then
      message "Incomplete: ($check_desc)" "$debug_out" -1
    else
      message "Missing: ($check_desc)" "$debug_out" -1
    fi
    return 0
  fi

  # Build optional --debug_out flag (omit entirely when debug_out is empty)
  local debug_flag=""
  [ -n "$debug_out" ] && debug_flag="--debug_out $debug_out"

  if [ ! -e "cfis" ]; then
    ln -sf $config_dir cfis
  fi

  command "update_runs_log_file.py" $dry_run

  # Run job script
  command "job_sp_canfar_v2.0.bash -c $config_dir -p $psf -r $retrieve --tile_det $tile_det --tile_mask $tile_mask -j $tile_job --n_smp $N_SMP --nsh_jobs $N_SMP $debug_flag" $dry_run 2>&1 | tee -a "$log_file"
}


if [ "$type" == "data" ]; then

    echo "Running on data"
    retrieve="vos"
    config_dir=$HOME/shapepipe/example/cfis
    export SP_DIR=$dir
    export SP_CONFIG=$config_dir

elif [ "$type" == "image_sims" ]; then

    echo "Running on image simulations"
    retrieve="symlink"
    config_dir=$HOME/shapepipe/example/cfis_image_sims
    # SP_DIR points to the run directory where input_tiles and input_exp live;
    # configs use $SP_DIR/input_* so those dirs stay outside SP_RUN and are
    # not found twice by ShapePipe's recursive glob scan.
    export SP_DIR=$dir
    export SP_CONFIG=$config_dir
    tile_det='sx'

else

    echo "Invalid input type $type" 

fi

echo "config_dir=$config_dir"


# Init message
if [ "$test_only" == "1" ]; then
  message "$(basename "$0") test mode, exiting." "$debug_out" 0
else
  message "$(basename "$0") processing mode, starting." "$debug_out" -1
fi


## Check options
message "checking options" "$debug_out" -1

if [ "$job" == "-1" ]; then
  message "No job indicated, use option -j" "$debug_out" 2
fi

if [ "$ID" == "-1" ]; then
  message "No image ID indicated, use option -e" "$debug_out" 3
fi

if [ "$psf" != "psfex" ] && [ "$psf" != "mccd" ]; then
  message "PSF (option -p) needs to be 'psfex' or 'mccd', not '$psf'" "$debug_out" 4
fi

if [ "$dry_run" != "0" ] && [ "$dry_run" != "1" ]; then
  message "dry_run must be 0 or 1, not $dry_run" "$debug_out" 8
fi


# Start script

source $HOME/shapepipe/scripts/sh/functions.sh

message "Starting $(basename "$0")" "$debug_out" -1
message "`date`" "$debug_out" -1
message "ID=$ID" "$debug_out" -1

if [ "$dry_run" == "1" ]; then
  message "running in dry run mode" "$debug_out" -1
fi

### PSF model letter: 'P' (psfex) or 'M' (mccd)
letter=${psf:0:1}
Letter=${letter^}

cd $dir

# Derive tile path components from ID (e.g. "000.227" -> IDra="000")
IDra=${ID%%.*}
work_dir="$dir/tiles/$IDra/$ID"
log_file="$work_dir/job_sp_canfar_v2.0.log"

# Create tile work directory
[ ! -d "$work_dir" ] && command "mkdir -p $work_dir" $dry_run
cd $work_dir
echo "$0 $@" > "$log_file"

# Write ID to first input
# Image sims use dash format (e.g. 233-293); real data uses dot format (233.293)
# which ShapePipe's in2out_pattern converts to dashes for output naming only,
# not for input file lookup — so write the format that matches the actual files.
if [ ! -e tile_numbers.txt ]; then
  if [ "$type" == "image_sims" ]; then
    echo ${ID//./-} > tile_numbers.txt
  else
    echo $ID > tile_numbers.txt
  fi
fi

# Output directory
if [ ! -d "output" ]; then
  command "mkdir output" $dry_run
fi


echo -n "pwd: "; pwd


# Avoid Qt error with setools
export DISPLAY=:1.0


# Run job — on scratch if available, otherwise in-place
# Job 1: download tile images (config_tile_Git_vos.ini)

if [ -n "$scratch" ]; then
  scratch_work="$scratch/tiles/$IDra/$ID"
  message "Copying work dir to scratch: $scratch_work" "$debug_out" -1
  [ ! -d "$scratch/tiles/$IDra" ] && command "mkdir -p $scratch/tiles/$IDra" $dry_run
  command "cp -aR . $scratch_work" $dry_run
  command "cd $scratch_work" $dry_run
fi

IDdec=${ID##*.}

(( do_job = job & 1 ))
if [[ $do_job != 0 ]]; then
  # Job 1: download tile images and weights
  if [ "$type" == "image_sims" ]; then
    n_exp=2
  else
    n_exp=4
  fi
  run_tile_job 1 "Git" "get_images_runner:${n_exp}"
fi

(( do_job = job & 2 ))
if [[ $do_job != 0 ]]; then
  if [ "$type" == "image_sims" ]; then
    # Image sims weights are already uncompressed; fake the Uz output directory
    # so downstream jobs can find the weight via last:uncompress_fits_runner.
    weight_src="$dir/input_tiles/CFIS_simu_weight-${ID//./-}.fits"
    if [ "$check" == "1" ]; then
      uz_run_dir=$(ls -dt "$work_dir/output/run_sp_tile_Uz"* 2>/dev/null | head -1)
      if [ -n "$uz_run_dir" ] && [ -e "$uz_run_dir/uncompress_fits_runner/output/$(basename $weight_src)" ]; then
        message "Complete: Uz $(basename $weight_src)" "$debug_out" -1
      else
        message "Missing: Uz $(basename $weight_src)" "$debug_out" -1
      fi
    else
      uz_out="$work_dir/output/run_sp_tile_Uz$(date +_%Y-%m-%d_%H-%M-%S)/uncompress_fits_runner/output"
      command "mkdir -p $uz_out" $dry_run
      if [ -e "$weight_src" ] && [ ! -e "$uz_out/$(basename $weight_src)" ]; then
        command "ln -sf $weight_src $uz_out/$(basename $weight_src)" $dry_run
      fi
    fi
  else
    # Job 2: uncompress tile weights
    run_tile_job 2 "Uz" "uncompress_fits_runner:1"
  fi
fi

(( do_job = job & 4 ))
if [[ $do_job != 0 ]]; then
  # Job 4: find exposures
  run_tile_job 4 "Fe" "find_exposures_runner:1"
fi

(( do_job = job & 8 ))
if [[ $do_job != 0 ]]; then
  # Job 8: retrieve exposure images
  if [ "$type" == "image_sims" ]; then
    n_exp=3
  else
    n_exp=6
  fi
  run_exp_job 8 "Gie" "get_images_runner:${n_exp}"
fi

(( do_job = job & 16 ))
if [[ $do_job != 0 ]]; then
  # Job 16: split exposures, get WCS headers
  run_exp_job 16 "Sp" "split_exp_runner:121"
fi

(( do_job = job & 32 ))
if [[ $do_job != 0 ]]; then
  # Job 32: mask exposures
  run_exp_job 32 "Ma" "mask_runner:40"
fi

(( do_job = job & 64 ))
if [[ $do_job != 0 ]]; then
  # Job 64: PSF model
  # For image_sims: fake PSF runs as part of job 512 (requires sexcat from job 256)
  # For data: run full exposure-level PSF modelling pipeline
  if [ "$type" == "image_sims" ]; then
    message "Job 64 (fake PSF) is handled as part of job 512 for image_sims — skipping." "$debug_out" -1
  elif [ "$psf" == "psfex" ]; then
    run_exp_job 64 "SxSePsf${Letter}i" "sextractor_runner:80 psfex_runner:80 psfex_interp_runner:40::warn setools_runner:80:rand_split"
  else
    message "MCCD not implemented yet for v2.0" "$debug_out" 10
  fi
fi

(( do_job = job & 128 ))
if [[ $do_job != 0 ]]; then
  # Job 128: merge exposure WCS headers into tile-level sqlite log
  run_tile_job 128 "Mh_exp" "merge_headers_runner:1"
fi

(( do_job = job & 256 ))
if [[ $do_job != 0 ]]; then
  # Job 256: object selection on tiles
  if [ "$tile_det" == "uc" ]; then
    run_tile_job 256 "Gic Uc" "get_images_runner:2 read_ext_sexcat_runner:1"
  else
    n_exp=2
    run_tile_job 256 "Sx" "sextractor_runner:$n_exp"
  fi
fi

(( do_job = job & 512 ))
if [[ $do_job != 0 ]]; then
  # Job 512: process tiles ([PSF interp,] vignets)
  # For image_sims: fake PSF runs first (requires sexcat from job 256), then vignets
  if [ "$type" == "data" ]; then
      run_tile_job 512 "${Letter}iViVi ${Letter}iViVi ${Letter}iViVi" "psfex_interp_runner:1 vignetmaker_runner_run_1:1 vignetmaker_runner_run_2:4"
  else
      run_tile_job 64 "fpsf" "fake_psf_runner:1"
      run_tile_job 512 "ViVi VViVi" "vignetmaker_runner_run_1:1 vignetmaker_runner_run_2:4"
  fi
fi

(( do_job = job & 1024 ))
if [[ $do_job != 0 ]]; then
  # Job 1024: shape measurement
  run_tile_job 1024 "Ng" "ngmix_runner:1"
fi

(( do_job = job & 2048 ))
if [[ $do_job != 0 ]]; then
  # Job 2048: merge catalogues
  run_tile_job 2048 "Mc_${psf}" "make_cat_runner:1"
fi


if [ -n "$scratch" ]; then
  message "Syncing output from scratch back to permanent dir" "$debug_out" -1
  command "rsync -a output/ $work_dir/output/" $dry_run
  command "cd $work_dir" $dry_run
  command "update_runs_log_file.py" $dry_run
  command "rm -rf $scratch_work" $dry_run
  command "cd $dir" $dry_run
fi

cd $dir

message "End $(basename "$0")" "$debug_out" -1
