#!/usr/bin/env bash

# Name: job_sp_multiband.sh
# Description: Process UNIONS tiles to
#              create multi-band catalogue
# Author: Martin Kilbinger, <martin.kilbinger@cea.fr>
# Date: 13/11/2020
# Package: ShapePipe

# Functions

## Run command. Stop script on error occurs, upload sp log files before stopping script.
command_sp() {
  cmd=$1

  $cmd
  res=$?

  RED='\033[0;31m'
  GREEN='\033[0;32m'
  NC='\033[0m' # No Color

  if [ $res != 0 ]; then
    echo -e "${RED}Error occured, exiting, '$cmd' returned $res${NC}"
    exit $res
  else
    echo -e "${GREEN}success, '$cmd' returned $res${NC}"
  fi
}

# Command line arguments

## Default values
job=31
survey='unions'
retrieve='vos'
do_env=0
match=1

## Help string
usage="Usage: $(basename "$0") [OPTIONS] TILE_ID_1 [TILE_ID_2 [...]]
\n\nOptions:\n
   -h\tthis message\n
   -e\tset environment and exit (run as '. $(basename "$0")'\n
   -j, --job JOB\tRunning JOB, bit-coded\n
   \t  1: retrieve images (online if method=vos)\n
   \t  2: prepare images (offline)\n
   \t  4: mask (online)\n
   \t  8: processing of stars on exposures (offline)\n
   \t 16: detection and matching (offline)\n
   \t 32: shapes and morphology (offline)\n
   \t 64: paste catalogues (offline)\n
   -s, --survey NAME\n
   \t survey name, one in ['unions'|'ps3pi_cfis']\n
   -r, --retrieve METHOD\n
   \tmethod to retrieve images, one in 'vos|symlink', default='$retrieve'\n
   --no_match\n
   \tdo not match with spectroscopic catalog\n
   TILE_ID_i\n
   \ttile ID, e.g. 282.247\n"

if [ -z $1 ]; then
        echo -ne $usage
        exit 1
fi

## Parse command line
ID=()
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
    -s|--survey)
      survey="$2"
      shift
      ;;
    -r|--retrieve)
      retrieve="$2"
      shift
      ;;
    --no_match)
      match=0
      ;;
    *)
      ID+=("$1")
      ;;
  esac
  shift
done

if [ $do_env == 1 ]; then
   echo "environment set, exiting"
   return
   exit 0
fi

if [ $match == 0 ]; then
  match_str="_nomatch"
else
  match_str=""
fi

n_ID=${#ID[@]}
echo "Processing $n_ID tile(s)"

# Path variables
export SP_HOME=$HOME/astro/repositories/github/shapepipe
export SP_CONFIG=$SP_HOME/example/$survey
export SP_RUN=.

# Create output path
mkdir -p output

# Write tile numbers to ASCII input file
rm -rf tile_numbers.txt
for TILE in ${ID[@]}; do
   echo $TILE >> tile_numbers.txt
done

# Processing: run pipeline

## Retrieve images
(( do_job= $job & 1 ))
if [[ $do_job != 0 ]]; then

  ### Retrieve tiles
  command_sp "shapepipe_run -c $SP_CONFIG/config_get_tiles_$retrieve.ini"

  ### Find exposures that were used to create tile (stack)
  command_sp "shapepipe_run -c $SP_CONFIG/config_tile_Fe.ini"

  ### Retrieve exposures
  command_sp "shapepipe_run -c $SP_CONFIG/config_exp_Gi_$retrieve.ini"

fi

## Prepare images
(( do_job= $job & 2 ))
if [[ $do_job != 0 ]]; then

  ### Convert Pan-STARRS image names, necessary for UNIONS
  ### but not PS3PI tiles
  if [ "$survey" == "unions" ]; then
    command_sp "$SP_HOME/scripts/python/ps_convert_file_names.py"
  fi

  ### Unzip FITS files and/or remove unused HDU:
  ### CFIS weights; UNIONS images, weights, flags
  command_sp "shapepipe_run -c $SP_CONFIG/config_unfz.ini"

  ### Split images into single-HDU files
  command_sp "shapepipe_run -c $SP_CONFIG/config_exp_Sp.ini"

  ### Merge FITS headers for WCS information
  command_sp "shapepipe_run -c $SP_CONFIG/config_exp_Mh.ini"

fi

## Mask
(( do_job= $job & 4 ))
if [[ $do_job != 0 ]]; then

  # Note: Before activating shapepipe, the following modules need to be loaded:
  #  intelpython/3
  #  openmpi/4.0.5

  ### Create flags for CFIS r-band images: add star, halo, and Messier
  ### object masks.
  command_sp "shapepipe_run -c $SP_CONFIG/config_tile_mask_r.ini"

  ### Mask r-band exposures
  command_sp "shapepipe_run -c $SP_CONFIG/config_exp_Ma.ini"

fi

## processing of stars on exposures (detection, and matching, PSF model)
(( do_job= $job & 8 ))
if [[ $do_job != 0 ]]; then

  ### Detect star candidates
  command_sp "shapepipe_run -c $SP_CONFIG/config_exp_Sx.ini"

  ### Select stars
  command_sp "shapepipe_run -c $SP_CONFIG/config_exp_Se.ini"

  ### Create  PSF model
  command_sp "shapepipe_run -c $SP_CONFIG/config_exp_Psm.ini"

  ### Bad hacks to get PSF input dir
  input_psfex=`find . -name star_selection-*.psf | head -n 1`
  command_sp "ln -sf `dirname $input_psfex` input_psfex"
  input_split_exp=`find output -name flag-*.fits | head -n 1`
  command_sp "ln -sf `dirname $input_split_exp` input_split_exp"
  input_sextractor=`find . -name sexcat_sexcat-*.fits | head -n 1`
  command_sp "ln -sf `dirname $input_sextractor` input_sextractor"

fi

## Detectino and matching
(( do_job= $job & 16 ))
if [[ $do_job != 0 ]]; then

  ### Detect objects on r-band images, measure properties
  ### on other bands
  command_sp "shapepipe_run -c $SP_CONFIG/config_tile_detect_r_me.ini"
  command_sp "shapepipe_run -c $SP_CONFIG/config_tile_detect_r.ini"
  command_sp "shapepipe_run -c $SP_CONFIG/config_tile_detect_i.ini"
  command_sp "shapepipe_run -c $SP_CONFIG/config_tile_detect_g.ini"
  command_sp "shapepipe_run -c $SP_CONFIG/config_tile_detect_z.ini"
  command_sp "shapepipe_run -c $SP_CONFIG/config_tile_detect_y.ini"

  if [ "$survey" == "unions" ]; then
    command_sp "shapepipe_run -c $SP_CONFIG/config_tile_detect_u.ini"
    command_sp "shapepipe_run -c $SP_CONFIG/config_tile_match_ext_u.ini"
  fi

  ### Match with external spectroscopic catalogue
  command_sp "shapepipe_run -c $SP_CONFIG/config_tile_match_ext_r_me.ini"
  command_sp "shapepipe_run -c $SP_CONFIG/config_tile_match_ext_r.ini"
  command_sp "shapepipe_run -c $SP_CONFIG/config_tile_match_ext_g.ini"
  command_sp "shapepipe_run -c $SP_CONFIG/config_tile_match_ext_i.ini"
  command_sp "shapepipe_run -c $SP_CONFIG/config_tile_match_ext_z.ini"
  command_sp "shapepipe_run -c $SP_CONFIG/config_tile_match_ext_y.ini"

  ## TODO: match_ext_u for unions

  ### Vignets for weights
  command_sp "shapepipe_run -c $SP_CONFIG/config_tile_Viw.ini"

  ### Interpolate exposure PSFs to tile objects
  command_sp "shapepipe_run -c $SP_CONFIG/config_tile_Psi.ini"

  ### Vignets for exposures
  command_sp "shapepipe_run -c $SP_CONFIG/config_tile_Vix.ini"

fi

## Shapes and morphology
(( do_job= $job & 32 ))
if [[ $do_job != 0 ]]; then

  command_sp "shapepipe_run -c $SP_CONFIG/config_tile_Sh${match_str}.ini"

fi

## Paste catalogs
(( do_job= $job & 64 ))
if [[ $do_job != 0 ]]; then

  command_sp "shapepipe_run -c $SP_CONFIG/config_tile_paste_cat_morph.ini"

fi

