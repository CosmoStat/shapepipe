#!/usr/bin/env bash

# Name: job_sp.bash
# Description: General script to process one or more tiles
#              with all contributing exposures.
#              This works as job submission script for
#              the canfar batch system.
#              called in interactive mode on a virtual
#              machine.
# Author: Martin Kilbinger <martin.kilbinger@cea.fr>
# Date: v1.0 11/2020
#       v1.1 01/2021

# MKDEBUG TODO:
# Option to change SMP_BATCH_SIZE, not for MPI


# VM home, required for canfar run.
## On other machines set to $HOME
export VM_HOME=/home/ubuntu
if [ ! -d "$VM_HOME" ]; then
    export VM_HOME=$HOME
fi

# Command line arguments
## Default values
job=255
#config_dir='vos:cfis/cosmostat/kilbinger/cfis'
config_dir=$VM_HOME/shapepipe/example/cfis
psf='mccd'
retrieve='vos'
star_cat_for_mask='onthefly'
results='cosmostat/kilbinger/results_v1'
n_smp=-1
nsh_step=-1
nsh_max=-1
nsh_jobs=8

## Help string
usage="Usage: $(basename "$0") [OPTIONS] TILE_ID_1 [TILE_ID_2 [...]]
\n\nOptions:\n
   -h\tthis message\n
   -j, --job JOB\tRunning JOB, bit-coded\n
   \t   1: retrieve images (online if method=vos)\n
   \t   2: prepare images (offline)\n
   \t   4: mask (online)\n
   \t   8: detection of galaxies on tiles; processing of stars on exposures (offline)\n
   \t  16: galaxy selection on tiles (offline)\n
   \t  32: shapes and morphology (offline)\n
   \t  64: paste catalogues (offline)\n
   \t 128: upload results (online)\n
   -c, --config_dir DIR\n
   \t config file directory, default='$config_dir'\n
   -p, --psf MODEL\n
    \tPSF model, one in ['psfex'|'mccd'], default='$psf'\n
   -r, --retrieve METHOD\n
   \tmethod to retrieve images, allowed are 'vos', 'symlink', default='$retrieve'\n
   -s, --star_cat_for_mask\n
   \tcatalogue for masking bright stars, allowed are 'onthefly', 'save',\n
   \tdefault is '${star_cat_for_mask}'\n
   -o, --output_dir\n
   \toutput (upload) directory on vos:cfis, default='$results'\n
   -n, --n_smp\n
   \tnumber of jobs (SMP mode only), default from original config files\n
   --nsh_step NSTEP\n
   --nsh_jobs NJOB\n
   \tnumber of shape measurement parallel jobs, default=$nsh_jobs\n
   \tnumber of objects per parallel shape module call, \n
   \tdefault: optimal number is computed\n
   --nsh_max NMAX\n
   \tmax number of objects per parallel shape module call, \n
   \tdefault: unlimited; has precedent over --nsh_step\n
   TILE_ID_i\n
   \ttile ID(s), e.g. 283.247 214.242\n
"

## Help if no arguments
if [ -z $1 ]; then
        echo -ne $usage
        exit 1
fi

## Parse command line
TILE_ARR=()
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
    -c|--config_dir)
      config_dir="$2"
      shift
      ;;
    -p|--psf)
      psf="$2"
      shift
      ;;
    -r|--retrieve)
      retrieve="$2"
      shift
      ;;
    -s|--star_cat_for_mask)
      star_cat_for_mask="$2"
      shift
      ;;
    -o|--output_dir)
      results="$2"
      shift
      ;;
    -n|--n_smp)
      n_smp="$2"
      shift
      ;;
    --nsh_max)
      nsh_max="$2"
      shift
      ;;
    --nsh_step)
      nsh_step="$2"
      shift
      ;;
    --nsh_jobs)
      nsh_jobs="$2"
      shift
      ;;
    *)
      TILE_ARR+=("$1")
      ;;
  esac
  shift
done

## Check options
if [ "$psf" != "psfex" ] && [ "$psf" != "mccd" ]; then
  echo "PSF (option -p) needs to be 'psfex' or 'mccd'"
  exit 2
fi

if [ "$star_cat_for_mask" != "onthefly" ] && [ "$star_cat_for_mask" != "save" ]; then
  echo "Star cat for mask (option -s) needs to be 'onthefly' or 'save'"
  exit 4
fi

if [ "$retrieve" != "vos" ] && [ "$retrieve" != "symlink" ]; then
  echo "method to retrieve images (option -r) needs to be 'vos' or 'symlink'"
  exit 5
fi

n_tile=${#TILE_ARR[@]}
if [ "$n_tile" == "0" ]; then
  echo "No tile ID given"
  exit 3
fi

if [ $nsh_max != -1 ]; then
  nsh_step=$nsh_max
fi

# For tar archives. Should be unique to each job
export ID=`echo ${TILE_ARR[@]} | tr ' ' '_'`

## Paths

# SExtractor library bug work-around
export PATH="$PATH:$VM_HOME/bin"

## Path variables used in shapepipe config files

# Run path and location of input image directories
export SP_RUN=`pwd`

# Config file path
export SP_CONFIG=$SP_RUN/cfis
export SP_CONFIG_MOD=$SP_RUN/cfis_mod

## Other variables

# Input tile numbers ASCII file
export TILE_NUMBERS_PATH=tile_numbers.txt

# Output
OUTPUT=$SP_RUN/output

# For tar archives
output_rel=`realpath --relative-to=. $OUTPUT`

# Stop on error, default=1
STOP=1

# Verbose mode (1: verbose, 0: quiet)
VERBOSE=1

# VCP options
export CERTFILE=$VM_HOME/.ssl/cadcproxy.pem
export VCP="vcp --certfile=$CERTFILE"


## Functions

# Print string, executes command, and prints return value.
function command () {
   cmd=$1
   str=$2

   RED='\033[0;31m'
   GREEN='\033[0;32m'
   NC='\033[0m' # No Color
   # Color escape characters show up in log files
   #RED=''
   #GREEN=''
   #NC=''


   if [ $# == 2 ]; then
      if [ $VERBOSE == 1 ]; then
           echo "$str: running '$cmd'"
      fi
      $cmd
   else
      if [ $VERBOSE == 1 ]; then
         echo "$str: running '$cmd $4 \"$5 $6\"'"
      fi
      $cmd $4 "$5 $6"
   fi	
   res=$?

   if [ $VERBOSE == 1 ]; then
      if [ $res == 0 ]; then
         echo -e "${GREEN}success, return value = $res${NC}"
      else
         echo -e "${RED}error, return value = $res${NC}"
         if [ $STOP == 1 ]; then
            echo "${RED}exiting 'canfar_sp.bash', error in command '$cmd'${NC}"
            exit $res
         else
            echo "${RED}continuing 'canfar_sp.bash', error in command '$cmd'${NC}"
         fi
      fi
   fi
}

# Run shapepipe command. If error occurs, upload sp log files before stopping script.
function command_sp() {
   local cmd=$1
   local str=$2

   command "$1" "$2"
}

# Set up config file and call shapepipe_run
function command_cfg_shapepipe() {
    local config_name=$1
    local str=$2
    local _n_smp=$3 

    config_upd=$(set_config_n_smp $config_name $_n_smp)
    local cmd="shapepipe_run -c $config_upd"
    command_sp "$cmd" "$str"
}

# Tar and upload files to vos
function upload() {
   base=$1
   shift
   ID=$1
   shift
   verbose=$1
   shift
   upl=("$@")

   echo "Counting upload files"
   n_upl=(`ls -l ${upl[@]} | wc`)
   if [ $n_upl == 0 ]; then
      if [ $STOP == 1 ]; then
         echo "Exiting script, no file found for '$base' tar ball"
         exit 3
      fi
   fi
   tar czf ${base}_${ID}.tgz ${upl[@]}
   command "$VCP ${base}_${ID}.tgz vos:cfis/$results" "Upload tar ball"
}

# Upload log files
function upload_logs() {
   id=$1
   verbose=$2

   upl="$output_rel/*/*/logs $output_rel/*/logs"
   upload "logs" "$id" "$verbose" "${upl[@]}"
}

# Print script variables
function print_env() {
   echo "*** Environment ***"
   echo "Data:"
   echo " TILE_ARR=${TILE_ARR[@]}"
   echo "Paths:"
   echo " VM_HOME=$VM_HOME"
   echo " SP_RUN=$SP_RUN"
   echo " TILE_NUMBERS_PATH=$TILE_NUMBERS_PATH"
   echo " OUTPUT=$OUTPUT"
   echo " SP_CONFIG=$SP_CONFIG"
   echo "Other variables:"
   echo " VCP=$VCP"
   echo " CERTFILE=$CERTFILE"
   echo "***"
}

function set_config_n_smp() {
  local config_name=$1
  local _n_smp=$2

  local config_orig="$SP_CONFIG/$config_name"

  if [[ $_n_smp != -1 ]]; then
    # Update SMP batch size
    local config_upd="$SP_CONFIG_MOD/$config_name"
    update_config $config_orig $config_upd "SMP_BATCH_SIZE" $_n_smp
  else
    # Keep original config file
    local config_upd=$config_orig
  fi

  # Set "return" value (stdout)
  echo "$config_upd"
}

# Update config file
function update_config() {
  local config_orig=$1
  local config_upd=$2
  local key=$3
  local val_upd=$4

  cat $config_orig \
    | perl -ane 's/'$key'\s+=.+/'$key' = '$val_upd'/; print' > $config_upd
}

### Start ###

echo "Start"

echo "Processing $n_tile tile(s)"

# Create input and output directories
mkdir -p $SP_RUN
cd $SP_RUN
mkdir -p $OUTPUT
mkdir -p $SP_CONFIG_MOD

# Processing

## Retrieve config files and images (online if retrieve=vos)
## Retrieve and save star catalogues for masking (if star_cat_for_mask=save)
(( do_job= $job & 1 ))
if [[ $do_job != 0 ]]; then

  # Write tile numbers to ASCII input file
  rm -rf $TILE_NUMBERS_PATH
  for TILE in ${TILE_ARR[@]}; do
    echo $TILE >> $TILE_NUMBERS_PATH
  done

  ### Retrieve config files
  if [[ $config_dir == *"vos:"* ]]; then
    command_sp "$VCP $config_dir ." "Retrieve shapepipe config files"
  else
    if [[ ! -L cfis ]]; then
      command_sp "ln -s $config_dir cfis" "Retrieve shapepipe config files"
    fi
  fi

  ### Retrieve files
  command_sp "shapepipe_run -c $SP_CONFIG/config_GitFeGie_$retrieve.ini" "Retrieve images"

  ### Retrieve and save star catalogues for masking
  if [ "$star_cat_for_mask" == "save" ]; then
    #### For tiles
    mkdir $SP_RUN/star_cat_tiles
    command_sp "create_star_cat $SP_RUN/output/run_sp_GitFeGie_*/get_images_runner_run_1/output $SP_RUN/star_cat_tiles" "Save star cats for masking (tile)"

    #### For single-exposures
    mkdir $SP_RUN/star_cat_exp
    command_sp "create_star_cat $SP_RUN/output/run_sp_GitFeGie_*/get_images_runner_run_2/output $SP_RUN/star_cat_exp exp" "Save star cats for masking (exp)"
  fi

fi

## Prepare images (offline)
(( do_job= $job & 2 ))
if [[ $do_job != 0 ]]; then

  ### Uncompress tile weights
  command_cfg_shapepipe "config_tile_Uz.ini" "Run shapepipe (uncompress tile weights)" $n_smp

  ### Split images into single-HDU files, merge headers for WCS info
  command_cfg_shapepipe "config_exp_SpMh.ini" "Run shapepipe (split images, merge headers)" $n_smp

fi

## Mask tiles and exposures: add star, halo, and Messier object masks (online if "star_cat_for_mask" is "onthefly")
(( do_job= $job & 4 ))
if [[ $do_job != 0 ]]; then

  ### Mask tiles and exposures
  command_cfg_shapepipe "config_MaMa_$star_cat_for_mask.ini" "Run shapepipe (mask)" $n_smp

fi


## Remaining exposure processing (offline)
(( do_job= $job & 8 ))
if [[ $do_job != 0 ]]; then

  ### Star detection, selection, PSF model. setools can exit with an error for CCD with insufficient stars,
  ### the script should continue
  STOP=0
  command_sp "shapepipe_run -c $SP_CONFIG/config_tile_Sx_exp_${psf}.ini" "Run shapepipe (tile detection, exp $psf)"
  STOP=1

fi

## Process tiles up to shape measurement
(( do_job= $job & 16 ))
if [[ $do_job != 0 ]]; then

  ### PSF model letter: 'P' (psfex) or 'M' (mccd)
  letter=${psf:0:1}
  Letter=${letter^}
  command_sp "shapepipe_run -c $SP_CONFIG/config_tile_${Letter}iViSmVi.ini" "Run shapepipe (tile PsfInterp=$Letter}: up to ngmix+galsim)"

fi

## Shape measurement (offline)
(( do_job= $job & 32 ))
if [[ $do_job != 0 ]]; then

  ### Prepare config files
  n_min=0
  if [[ $nsh_step == -1 ]]; then
    n_obj=`get_number_objects.py`
    nsh_step=`echo "$(($n_obj/$nsh_jobs))"`
  fi

  n_max=$((nsh_step - 1))
  for k in $(seq 1 $nsh_jobs); do
    cat $SP_CONFIG/config_tile_Ng_template.ini | \
      perl -ane \
        's/(ID_OBJ_MIN =) X/$1 '$n_min'/; s/(ID_OBJ_MAX =) X/$1 '$n_max'/; s/NgXu/Ng'$k'u/; s/X_interp/'$psf'_interp/g; print' \
        > $SP_CONFIG_MOD/config_tile_Ng${k}u.ini
    n_min=$((n_min + nsh_step))
    if [ "$k" == $((nsh_jobs - 1)) ] && [ $nsh_max == -1 ]; then
      n_max=-1
    else
      n_max=$((n_min + nsh_step - 1))
    fi
  done

  ### Shapes, run $nsh_jobs parallel processes
  VERBOSE=0
  for k in $(seq 1 $nsh_jobs); do
      command_sp "shapepipe_run -c $SP_CONFIG_MOD/config_tile_Ng${k}u.ini" "Run shapepipe (tile: ngmix+galsim $k)" &
  done
  wait
  VERBOSE=1

fi

## Create final catalogues (offline)
(( do_job= $job & 64 ))
if [[ $do_job != 0 ]]; then

  cat $SP_CONFIG/config_merge_sep_cats_template.ini | \
    perl -ane \
      's/(N_SPLIT_MAX =) X/$1 '$nsh_jobs'/; print' \
      > $SP_CONFIG_MOD/config_merge_sep_cats.ini
 
  ### Merge separated shapes catalogues
  command_sp "shapepipe_run -c $SP_CONFIG_MOD/config_merge_sep_cats.ini" "Run shapepipe (tile: merge sep cats)" "$VERBOSE" "$ID"

  ### Merge all relevant information into final catalogue
  command_sp "shapepipe_run -c $SP_CONFIG/config_make_cat_$psf.ini" "Run shapepipe (tile: create final cat $psf)" "$VERBOSE" "$ID"

fi

## Upload results (online)
(( do_job= $job & 128 ))
if [[ $do_job != 0 ]]; then

  ### module and pipeline log files
  upload_logs "$ID" "$VERBOSE"

  ### Final shape catalog
  ### pipeline_flags are the tile masks, for random cats
  ### SETools masks (selection), stats and plots
  ### ${psf}_interp_exp for diagnostics, validation with leakage,
  ### validation with residuals, rho stats

  NAMES=(
    "final_cat"
    "pipeline_flag"
    "setools_mask"
    "setools_stat"
    "setools_plot"
  )
  DIRS=(
    "*/make_cat_runner/output"
    "*/mask_runner_run_1/output"
    "*/setools_runner/output/mask"
    "*/setools_runner/output/stat"
    "*/setools_runner/output/plot"
  )
  PATTERNS=(
    "final_cat-*"
    "pipeline_flag-???-???*"
    "*"
    "*"
    "*"
  )

  # PSF validation
  pattern="validation_psf-*"
  if [ "$psf" == "psfex" ]; then
    name="psfex_interp_exp"
    dir="*/psfex_interp_runner/output"
  else
    name="mccd_fit_val_runner"
    dir="*/mccd_fit_val_runner/output"
  fi
  upl=$output_rel/$dir/$pattern
  upload "$name" "$ID" "$VERBOSE" "${upl[@]}"

  for n in "${!NAMES[@]}"; do
      name=${NAMES[$n]}
      dir=${DIRS[$n]}
      pattern=${PATTERNS[$n]}
      upl=$output_rel/$dir/$pattern
      upload "$name" "$ID" "$VERBOSE" "${upl[@]}"
  done

fi
