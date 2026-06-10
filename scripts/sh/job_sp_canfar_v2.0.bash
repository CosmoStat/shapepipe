#!/bin/bash

# Name: job_sp_canfar.bash
# Description: General script to process one or more tiles
#              with all contributing exposures.
#              This works as job submission script for
#              the canfar batch system.
#              called in interactive mode on a virtual
#              machine.
# Author: Martin Kilbinger <martin.kilbinger@cea.fr>

# Shared job-list description
source $HOME/shapepipe/scripts/sh/job_list_help.bash

# Command line arguments
## Default values
job=255
config_dir=$HOME/shapepipe/example/cfis
psf='psfex'
retrieve='vos'
star_cat_for_mask='onthefly'
tile_det='sx'
tile_mask=0
exclusive=''
n_smp=-1
nsh_jobs=8
debug_out=""

pat="--- "

## Help string
usage="Usage: $(basename "$0") [OPTIONS] [TILE_ID]
\n\nOptions:\n
   -h\tthis message\n
   -j, --job JOB\tRunning JOB, bit-coded\n
${JOB_LIST_HELP}
   -c, --config_dir DIR\n
   \t config file directory, default='$config_dir'\n
   -p, --psf MODEL\n
    \tPSF model, one in ['psfex'|'mccd'], default='$psf'\n
   -r, --retrieve METHOD\n
   \tmethod to retrieve images, allowed are 'vos', 'symlink', default='$retrieve'\n
   -s, --star_cat_for_mask\n
   \tcatalogue for masking bright stars, allowed are 'onthefly', 'save',\n
   \tdefault is '${star_cat_for_mask}'\n
   --tile_det DET\n
   \ttile object detection mode, one in ['sx'|'uc'], default='${tile_det}'\n
   --tile_mask MASK\n
   \ttile masking, default='${tile_mask}'\n
   \t  sx: run SExtractor on tile image\n
   \t  uc: use UNIONS catalogue (Stephen Gwyn)\n
   -e, --exclusive ID\n
   \texclusive input filer number string ID (default: None)\n
   -n, --n_smp N_SMP\n
   \tnumber of jobs (SMP mode only), default from original config files\n
   --nsh_jobs NJOB\n
   \tnumber of objects per parallel shape module call, \n
   \tdefault: optimal number is computed\n
   --debug_out PATH\n
   \tdebug output file PATH, default not used\n
   TILE_ID_i\n
   \ttile ID(s), e.g. 283.247 214.242, only with '-j 1'\n
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
    --tile_det)
      tile_det="$2"
      shift
      ;;
    --tile_mask)
      tile_mask="$2"
      shift
      ;;
    -e|--exclusive)
      exclusive="$2"
      shift
      ;;
    -n|--n_smp)
      n_smp="$2"
      shift
      ;;
    --nsh_jobs)
      nsh_jobs="$2"
      shift
      ;;
    --debug_out)
      debug_out="$2"
      shift
      ;;
  esac
  shift
done

## Check options
if [ "$psf" != "psfex" ] && [ "$psf" != "mccd" ] && [ "$psf" != "psf" ]; then
  echo "PSF (option -p) needs to be 'psfex', 'mccd', or 'psf' (image sims)"
  exit 2
fi

if [ "$star_cat_for_mask" != "onthefly" ] && [ "$star_cat_for_mask" != "save" ]; then
  echo "Star cat for mask (option -s) needs to be 'onthefly' or 'save'"
  exit 4
fi

if [ "$retrieve" != "vos" ] && [ "$retrieve" != "symlink" ]; then
  echo "Invalid method to retrieve images $retrieve (option -r), needs to be 'vos' or 'symlink'"
  exit 5
fi

if [ "$tile_det" != "sx" ] && [ "$tile_det" != "uc" ]; then
  echo "tile detection mode (option --tile_det) needs to be 'sx' or 'uc'"
  exit 6
fi

if [ "$tile_mask" != 1 ] && [ "$tile_mask" != 0 ]; then
  echo "tile mask needs to be 0 or 1"
  exit 7
fi

if [ -n "$debug_out" ]; then
  echo $pat`date` >> $debug_out
  echo "${pat}Starting $(basename "$0")" >> $debug_out
fi

## Paths

## Path variables used in shapepipe config files

# Run path and location of input image directories
export SP_RUN=`pwd`

# Config file path — use value exported by run_job_sp_canfar_v2.0.bash if set,
# otherwise fall back to the cfis symlink in the run directory.
export SP_CONFIG=${SP_CONFIG:-$SP_RUN/cfis}

# Root directory for per-exposure work directories.
# Set SP_EXP in the environment to override; otherwise use SP_DIR (the run
# root, always exported by run_job_sp_canfar_v2.0.bash for both data and
# image_sims) so exp/ is always a sibling of tiles/ under the same root.
if [ -z "${SP_EXP}" ]; then
  export SP_EXP="$SP_DIR/exp"
  echo "Setting SP_EXP to $SP_EXP"
fi

## Other variables

# Override OMP_NUM_THREADS if the CANFAR provisioning template was not expanded
if [[ "${OMP_NUM_THREADS}" == *'${'* ]] || [[ "${OMP_NUM_THREADS}" == *'.'* ]]; then
  export OMP_NUM_THREADS=1
fi

# Output
OUTPUT=$SP_RUN/output

# Stop on error, default=1
STOP=1

# Verbose mode (1: verbose, 0: quiet)
VERBOSE=1


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


   if [ -n "$debug_out" ]; then
      echo "${pat}pwd = `pwd`" >> $debug_out
      echo "${pat}SP_RUN = $SP_RUN" >> $debug_out
      echo "${pat}SP_CONFIG = $SP_CONFIG" >> $debug_out
    fi

   if [ $VERBOSE == 1 ]; then
      echo "$str: running '$cmd'"
   fi

   if [ -n "$debug_out" ]; then
      echo "${pat}Running $cmd" >> $debug_out
   fi

      $cmd
      
   res=$?

   if [ -n "$debug_out" ]; then
       echo "${pat}exit code = $res" >> $debug_out
   fi

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

# Set up config file and call shapepipe_run.
# Batch size is passed via --batch_size flag; no config editing needed.
function command_cfg_shapepipe() {
    local config_name=$1
    local str=$2
    local _n_smp=$3
    local _exclusive=$4

    if [ "$exclusive" != "" ]; then
      exclusive_flag="-e $_exclusive"
    else
      exclusive_flag=""
    fi

    local batch_flag=""
    if [[ $_n_smp != -1 ]]; then
      batch_flag="--batch_size $_n_smp"
    fi

    local config="$SP_CONFIG/$config_name"
    local cmd="shapepipe_run -c $config $exclusive_flag $batch_flag"
    command "$cmd" "$str"
}

### Start ###

echo "Start processing"

# Create input and output directories
mkdir -p $SP_RUN
cd $SP_RUN
mkdir -p $OUTPUT

# Processing


## Retrieve tile images (online if retrieve=vos)
## Retrieve and save star catalogues for masking (if star_cat_for_mask=save)
(( do_job = $job & 1 ))
if [[ $do_job != 0 ]]; then

  ### Retrieve files
  command_cfg_shapepipe \
    "config_tile_Git_$retrieve.ini" \
     "Retrieve images" \
     -1 \
     $exclusive

  ### Retrieve and save star catalogues for masking
  if [ "$star_cat_for_mask" == "save" ]; then
    #### For tiles
    mkdir $SP_RUN/star_cat_tiles
    command \
      "create_star_cat $SP_RUN/output/run_sp_tile_Git_*/get_images_runner_run/output $SP_RUN/star_cat_tiles" \
      "Save star cats for masking (tile)"
  fi

fi

## Uncompress tile weights
(( do_job = $job & 2 ))
if [[ $do_job != 0 ]]; then

  command_cfg_shapepipe "config_tile_Uz.ini" "Run shapepipe (uncompress tile weights)" $n_smp $exclusive

fi

## Find exposures
(( do_job = $job & 4 ))
if [[ $do_job != 0 ]]; then

  command_cfg_shapepipe "config_tile_Fe.ini" "Run shapepipe (find exposures)" $n_smp $exclusive

fi

## Retrieve exposure images (online if retrieve=vos)
(( do_job = $job & 8 ))
if [[ $do_job != 0 ]]; then

  command_cfg_shapepipe \
    "config_exp_Gie_$retrieve.ini" \
    "Run shapepipe (get exposure images)" \
    $n_smp \
    $exclusive

fi

## Split images into single-HDU files, merge headers for WCS info
(( do_job = $job & 16 ))
if [[ $do_job != 0 ]]; then

  command_cfg_shapepipe \
    "config_exp_Sp.ini" \
    "Run shapepipe (split images, merge headers)" \
    $n_smp \
    $exclusive

fi

## Mask exposures
(( do_job = $job & 32 ))
if [[ $do_job != 0 ]]; then

  command_cfg_shapepipe \
    "config_exp_Ma_${star_cat_for_mask}.ini" \
    "Run shapepipe (mask exposures)" \
    $n_smp \
    $exclusive

fi

## Exposure processing (offline)
(( do_job = $job & 64 ))
if [[ $do_job != 0 ]]; then

  ### Star detection, selection, PSF model. setools can exit with an error for CCD with insufficient stars,
  ### the script should continue
  STOP=0
  command_cfg_shapepipe \
    "config_exp_${psf}.ini" \
    "Run shapepipe (exp $psf)" \
    $n_smp \
    $exclusive
  STOP=1

fi

## Merge single-exposure WCS headers into tile-level sqlite log.
(( do_job = $job & 128 ))
if [[ $do_job != 0 ]]; then

  command_cfg_shapepipe \
    "config_tile_Mh_exp.ini" \
    "Run shapepipe (merge exp headers)" \
    $n_smp \
    $exclusive

fi

## Select objects on tiles
(( do_job = $job & 256 ))
if [[ $do_job != 0 ]]; then

  if [ "$tile_det" == "uc" ]; then

    ### Download external catalogue from vos
    command_cfg_shapepipe \
      "config_tile_Git_cat_$retrieve.ini" \
      "Run shapepipe (download external tile catalogue)" \
      -1 \
      $exclusive

    ### Object selection from external catalogue
    command_cfg_shapepipe \
      "config_tile_Uc.ini" \
      "Run shapepipe (tile object selection, external catalogue)" \
      $n_smp \
      $exclusive

  else

    if [ "$tile_mask" == 0 ]; then

      ### Object detection on tiles with SExtractor
      command_cfg_shapepipe \
        "config_tile_Sx_nomask.ini" \
        "Run shapepipe (tile detection, SExtractor without input flags)" \
        $n_smp \
        $exclusive

    else

      ### Tile masking
      command_cfg_shapepipe \
        "config_tile_Ma_${star_cat_for_mask}.ini" \
        "Run shapepipe (tile detection, SExtractor with input flags)" \
        $n_smp \
        $exclusive

    ### Object detection on tiles with SExtractor
      command_cfg_shapepipe \
        "config_tile_Sx.ini" \
        "Run shapepipe (tile detection, SExtractor with input flags)" \
        $n_smp \
        $exclusive

    fi

  fi

fi

## Process tiles up to shape measurement: postage stamp creation
(( do_job = $job & 512 ))
if [[ $do_job != 0 ]]; then

  ### PSF model letter: 'P' (psfex) or 'M' (mccd)
  letter=${psf:0:1}
  Letter=${letter^}
  command_cfg_shapepipe \
    "config_tile_${Letter}iViVi_canfar_${tile_det}.ini" \
    "Run shapepipe (tile PsfInterp=${Letter}: postage stamp creation)" \
    $n_smp \
    $exclusive

fi

## Process tiles up to shape measurement: postage stamp creation
(( do_job = $job & 1024 ))
if [[ $do_job != 0 ]]; then

  command_cfg_shapepipe \
    "config_tile_Ng_batch_${psf}_${tile_det}.ini" \
    "Run shapepipe (tile): shape measurement" \
    $n_smp \
    $exclusive

fi

## Create final catalogues (offline)
(( do_job = $job & 2048 ))
if [[ $do_job != 0 ]]; then

  ### Merge all relevant information into final catalogue
  command_cfg_shapepipe \
    "config_tile_Mc_$psf.ini" \
    "Run shapepipe (tile: create final cat $psf)" \
    $n_smp \
    $exclusive

fi

if [ -n "$debug_out" ]; then
  echo "${pat}End $(basename "$0") ID=$exclusive success" >> $debug_out
fi
