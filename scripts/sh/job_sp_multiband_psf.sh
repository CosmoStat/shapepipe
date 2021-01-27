#!/usr/bin/env bash

# Name: job_sp_multiband.sh
# Description: Process UNIONS tiles to
#              create multi-band catalogue
# Author: Martin Kilbinger, <martin.kilbinger@cea.fr>
# Date: 13/11/2020
# Package: ShapePipe


usage="Usage: $(basename "$0") [OPTIONS] TILE_ID_1 [TILE_ID_2 [...]]
\n\nOptions:\n
   -h\tThis message\n
   -j, --job JOB\tRunning JOB, bit-coded\n
   \t  1: Retrieve images (online if method=vos)\n
   \t  2: Prepare images (offline)\n
   \t  4: Mask (online)\n
   \t  8: Remaining processing (offline)\n
   TILE_ID_i\n
   \ttile ID, e.g. 282.247\n"

if [ -z $1 ]; then
        echo -ne $usage
        exit 1
fi


# Command line arguments

## Default values
job=7

## Parse command line
ID=()
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
    *)
      ID+=("$1")
      ;;
  esac
  shift
done

n_ID=${#ID[@]}
echo "Processing $n_ID tile(s)"

# Path variables
export SP_HOME=$HOME/astro/repositories/github/shapepipe
export SP_CONFIG=$SP_HOME/example/unions
export SP_RUN=.

#export LD_LIBRARY_PATH=$HOME/.conda/envs/shapepipe/lib

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
  shapepipe_run -c $SP_CONFIG/config_get_tiles_symlink.ini

  ### Find exposures that were used to create tile (stack)
  shapepipe_run -c $SP_CONFIG/config_tile_Fe.ini

  ### Retrieve exposures
  shapepipe_run -c $SP_CONFIG/config_exp_Gi_symlink.ini

fi

## Prepare images
(( do_job= $job & 2 ))
if [[ $do_job != 0 ]]; then

  ### Convert Pan-STARRS image names, necessary for UNIONS
  ### but not PS3PI tiles
  $SP_HOME/scripts/python/ps_convert_file_names.py

  ### Unzip FITS files and/or remove unused HDU:
  ### CFIS weights; UNIONS images, weights, flags
  shapepipe_run -c $SP_CONFIG/config_unfz.ini

  ### Split images into single-HDU files
  shapepipe_run -c $SP_CONFIG/config_exp_Sp.ini

  ### Merge FITS headers for WCS information
  shapepipe_run -c $SP_CONFIG/config_exp_Mh.ini

fi

## Mask
(( do_job= $job & 4 ))
if [[ $do_job != 0 ]]; then

  # Note: Before activating shapepipe, the following modules need to be loaded:
  #  intelpython/3
  #  openmpi/4.0.5

  ### Create flags for CFIS r-band images: add star, halo, and Messier
  ### object masks.
  $CONDA_PREFIX/bin/mpiexec -np 4 shapepipe_run -c $SP_CONFIG/config_tile_mask_r.ini

  ### Mask r-band exposures
  $CONDA_PREFIX/bin/mpiexec -np 4 shapepipe_run -c $SP_CONFIG/config_exp_Ma.ini

fi

## Remaining processing
(( do_job= $job & 8 ))
if [[ $do_job != 0 ]]; then

  ### Detect star candidates
  shapepipe_run -c $SP_CONFIG/config_exp_Sx.ini

  ### Select stars
  shapepipe_run -c $SP_CONFIG/config_exp_Se.ini

  ### Create  PSF model
  shapepipe_run -c $SP_CONFIG/config_exp_Psm.ini

  ### Detect objects on r-band images, measure properties
  ### on u-, r-, i-, z-band images
  shapepipe_run -c $SP_CONFIG/config_tile_detect_r_me.ini
  shapepipe_run -c $SP_CONFIG/config_tile_detect_u.ini
  shapepipe_run -c $SP_CONFIG/config_tile_detect_i.ini
  shapepipe_run -c $SP_CONFIG/config_tile_detect_z.ini

  ### Vignets for weights
  shapepipe_run -c $SP_CONFIG/config_tile_Viw.ini

  ### Bad hack to get PSF input dir
  input_psfex=`find . -name star_selection-*.psf | head -n 1`
  ln -sf `dirname $input_psfex` input_psfex
  input_split_exp=`find output -name flag-*.fits | head -n 1`
  ln -sf `dirname $input_split_exp` input_split_exp
  input_sextractor=`find . -name sexcat_sexcat-*.fits | head -n 1`
  ln -sf `dirname $input_sextractor` input_sextractor

  ### Interpolate exposure PSFs to tile objects
  shapepipe_run -c $SP_CONFIG/config_tile_Psi.ini 

  ### Vignets for exposures
  shapepipe_run -c $SP_CONFIG/config_tile_Vix.ini

fi

(( do_job= $job & 16 ))
if [[ $do_job != 0 ]]; then

  ### Shapes and morphology
  shapepipe_run -c $SP_CONFIG/config_tile_Sh.ini

fi

(( do_job= $job & 32 ))
if [[ $do_job != 0 ]]; then

  # Merge catalogs
  shapepipe_run -c $SP_CONFIG/config_tile_paste_cat_morph.ini

fi

