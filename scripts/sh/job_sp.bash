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
# Package: shapepipe

# Command line arguments
## Default values
do_env=0
job=255
psf='mccd'
retrieve='vos'
nsh_step=3500
nsh_max=-1
nsh_jobs=8

## Help string
usage="Usage: $(basename "$0") [OPTIONS] TILE_ID_1 [TILE_ID_2 [...]]
\n\nOptions:\n
   -h\tthis message\n
   -e\tset environment and exit (run as '. $(basename "$0")'\n
   -j, --job JOB\tRunning JOB, bit-coded\n
   \t   1: retrieve images (online if method=vos)\n
   \t   2: prepare images (offline)\n
   \t   4: mask (online)\n
   \t   8: processing of stars on exposures (offline)\n
   \t  16: detection and vignets (offline)\n
   \t  32: shapes and morphology (offline)\n
   \t  64: paste catalogues (offline)\n
   \t 128: upload results (online)\n
   -p, --psf MODEL\n
    \tPSF model, one in ['psfex'|'mccd'], default='$psf'\n
   -r, --retrieve METHOD\n
   \tmethod to retrieve images, one in ['vos'|'symlink]', default='$retrieve'\n
   --nsh_jobs NJOB\n
   \tnumber of shape measurement parallel jobs, default=$nsh_jobs\n
   --nsh_step NSTEP\n
   \tnumber of objects per parallel shape module call, \n
   \t default: $nsh_step\n
   --nsh_max NMAX\n
   \tmax number of objects per parallel shape module call, \n
   \t default: unlimited; has precedent over --nsh_step\n
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
    -e)
      do_env=1
      ;;
    -j|--job)
      job="$2"
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

# VM home, required for canfar run.
# On other machines set to $HOME
export VM_HOME=/home/ubuntu
if [ ! -d "$VM_HOME" ]; then
    export VM_HOME=$HOME
fi

# SExtractor library bug work-around
export PATH="$PATH:$VM_HOME/bin"

# Results upload subdirectory on vos
RESULTS=results_$psf

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

# Stop on error
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

   #RED='\033[0;31m'
   #GREEN='\033[0;32m'
   #NC='\033[0m' # No Color
   # Color escape characters show up in log files
   RED=''
   GREEN=''
   NC=''


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

   #return $res
}

# Run shapepipe command. If error occurs, upload sp log files before stopping script.
command_sp() {
   cmd=$1
   str=$2

   #STOP=0
   command "$1" "$2"
   #res=$?
   #if [ $res != 0 ]; then
      #upload_logs $ID $VERBOSE
      #echo "exiting 'canfar_sp.bash', '$cmd' returned $res, log files for id=$ID uploaded"
      #exit $res
   #fi
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

   n_upl=(`ls -l ${upl[@]} | wc`)
   if [ $n_upl == 0 ]; then
      if [ $STOP == 1 ]; then
         echo "Exiting script, no file found for '$base' tar ball"
         exit 3
      fi
   fi
   tar czf ${base}_${ID}.tgz ${upl[@]}
   command "$VCP ${base}_${ID}.tgz vos:cfis/cosmostat/kilbinger/$RESULTS" "Upload log tar ball"
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
   echo " STOP=$STOP"
   echo " verbose=$VERBOSE"
   echo "***"
}


### Start ###

# Activate conda environment
if [ "$CONDA_DEFAULT_ENV" != "shapepipe" ]; then
  echo "Activate conda 'shapepipe' environment"
  source $VM_HOME/miniconda3/bin/activate shapepipe
fi

print_env

if [ $do_env == 1 ]; then
   echo "Exiting"
   return
fi

echo "Start"

echo "Processing $n_tile tile(s)"

# Create input and output directories
mkdir -p $SP_RUN
cd $SP_RUN
mkdir -p $OUTPUT

# Processing

## Retrieve config files and images (online if retrieve=vos)
(( do_job= $job & 1 ))
if [[ $do_job != 0 ]]; then

  # Write tile numbers to ASCII input file
  rm -rf $TILE_NUMBERS_PATH
  for TILE in ${TILE_ARR[@]}; do
    echo $TILE >> $TILE_NUMBERS_PATH
  done

  ### Download config files
  command_sp "$VCP vos:cfis/cosmostat/kilbinger/cfis ." "Get shapepipe config files"

  ### Get tiles
  command_sp "shapepipe_run -c $SP_CONFIG/config_get_tiles_$retrieve.ini" "Run shapepipe (get tiles)"

  ### Find exposures
  command_sp "shapepipe_run -c $SP_CONFIG/config_find_exp.ini" "Run shapepipe (find exposures)" 

  ### Retrieve exposures
  command_sp "shapepipe_run -c $SP_CONFIG/config_get_exp_$retrieve.ini" "Run shapepipe (get exposures)"

fi

## Prepare images (offline)
(( do_job= $job & 2 ))
if [[ $do_job != 0 ]]; then

  ### Uncompress tile weights
  command_sp "shapepipe_run -c $SP_CONFIG/config_unfz_w.ini" "Run shapepipe (uncompress tile weights)"

  ### Split images into single-HDU files, merge headers for WCS info
  command_sp "shapepipe_run -c $SP_CONFIG/config_exp_SpMh.ini" "Run shapepipe (split images, merge headers)"

fi

## Mask tiles and exposures: add star, halo, and Messier object masks (online)
(( do_job= $job & 4 ))
if [[ $do_job != 0 ]]; then

  ### Mask tiles
  command_sp "shapepipe_run -c $SP_CONFIG/config_tile_Ma.ini" "Run shapepipe (mask tiles)"

  ### Mask exposures
  command_sp "shapepipe_run -c $SP_CONFIG/config_exp_Ma.ini" "Run shapepipe (mask exposures)"

fi


## Remaining exposure processing (offline)
(( do_job= $job & 8 ))
if [[ $do_job != 0 ]]; then

  ### Star detection, selection, PSF model. setools can exit with an error for CCD with insufficient stars,
  ### the script should continue
  STOP=0
  command_sp "shapepipe_run -c $SP_CONFIG/config_exp_$psf.ini" "Run shapepipe (exp $psf)"
  STOP=1

  ### The following are very a bad hacks to get additional input file paths
  if [ "$psf" == "psfex" ]; then
    #input_psfex=`find . -name star_split_ratio_80-*.psf | head -n 1`
    input_psfex=`ls -1 ./output/*/psfex_runner/output/star_split_ratio_80*.psf | head -n 1`
    command_sp "ln -s `dirname $input_psfex` input_psfex" "Link psfex output"
  else
    input_psf_mccd=`find . -name "fitted_model*.npy" | head -n 1`
    command_sp "ln -s `dirname $input_psf_mccd` input_psf_mccd" "Link MCCD output"
  fi

  #input_split_exp=`find output -name flag-*.fits | head -n 1`
  input_split_exp=`ls -1 ./output/*/split_exp_runner/output/flag*.fits | head -n 1`
  command_sp "ln -s `dirname $input_split_exp` input_split_exp" "Link split_exp output"

  #input_sextractor=`find . -name sexcat_sexcat-*.fits | head -n 1`
  input_sextractor=`ls -1 ./output/*/sextractor_runner/output/sexcat_sexcat*.fits | head -n 1`
  command_sp "ln -s `dirname $input_sextractor` input_sextractor" "Link sextractor output"

fi

## Process tiles up to shape measurement
(( do_job= $job & 16 ))
if [[ $do_job != 0 ]]; then

  ### PSF model letter: 'P' (psfex) or 'M' (mccd)
  letter=${psf:0:1}
  Letter=${letter^}
  command_sp "shapepipe_run -c $SP_CONFIG/config_tile_Sx${Letter}iViSmVi.ini" "Run shapepipe (tile PsfInterp=$Letter}: up to ngmix+galsim)"

fi

## Shape measurement (offline)
(( do_job= $job & 32 ))
if [[ $do_job != 0 ]]; then

  ### Prepare config files
  mkdir -p $SP_CONFIG_MOD
  n_min=0
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

  ### ${psf}_interp_exp for diagnostics, validation with leakage,
  ### validation with residuals, rho stats
  ### SETools masks (selection), stats and plots
  ### pipeline_flags are the tile masks, for random cats
  ### Final shape catalog

  NAMES=(
    "final_cat"
    "pipeline_flag"
    "setools_mask"
    "setools_stat"
    "setools_plot"
  )
  DIRS=(
    "*/make_catalog_runner/output"
    "*/mask_runner/output"
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
