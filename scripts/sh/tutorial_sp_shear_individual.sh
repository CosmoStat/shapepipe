#!/usr/bin/env bash

# Name: tutorial_sp_shear_individual.sh
# Author: Martin Kilbinger <martin.kilbinger@cea.fr>
# Date: Dec 2020
# Descriptio: Runs a ShapePipe tutorial job
# Package: ShapePipe


# Default command-line variables
retrieve="vos"
SP_BASE=$HOME/astro/repositories/github/shapepipe

usage="Usage: $(basename "$0") [OPTIONS] TILE_ID_1 [TILE_ID_2 [...]]
\n\nOptions:\n
   -h\tThis message\n
   -r, --retrieve METHOD\n
   \tmethod to retrieve images, one in 'vos|symlink', default='$retrieve'\n
   -b, --base_dir_sp PATH\n
   \tShapePipe base directory, default='$SP_BASE'\n
   TILE_ID_i\n
   \ttile ID, e.g. 282.247\n"

if [ -z $1 ]; then
        echo -ne $usage
        exit 1
fi

# Command line arguments


# Parse command line
ID=()
while [ $# -gt 0 ]; do
  case "$1" in
    -h)
      echo -ne $usage
      exit 0
      ;;
    -r|--retrieve)
      retrieve="$2"
      shift
      ;;
    -b|--base_path_sp)
      SP_BASE="$2"
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
export SP_RUN=`pwd`
export SP_CONFIG=$SP_BASE/example/tutorial

# Create output path
mkdir -p output

# Pre-processing

## Select tile IDs to process
#cfis_field_select -i $SP_BASE/aux/CFIS/tiles_202007/tiles_all_order.txt --coord 213.68deg_54.79deg -t tile --input_format ID_only --out_name_only --out_ID_only -s -o tile_numbers -v

## Retrieve tiles
if [ "$retrieve" == "vos" ]; then

  shapepipe_run -c $SP_CONFIG/config_tile_Git.ini

elif [ "$retrieve" == "symlink" ]; then

  if [ ! -d "data" ]; then
    echo "Input directory 'data' not found"
    exit 2
  fi
  shapepipe_run -c $SP_CONFIG/config_tile_Git_symlink.ini

else

   echo "Invalid retrieve value '$retrieve'"

fi

## Uncompress tile weights
shapepipe_run -c $SP_CONFIG/config_tile_Uz.ini

## Find single exposures
shapepipe_run -c $SP_CONFIG/config_tile_Fe.ini

## Retrieve single exposures
shapepipe_run -c $SP_CONFIG/config_tile_Gie.ini

# Processing of single exposures

## Split into single-exposure single-HDU files and write FITS headers
shapepipe_run -c $SP_CONFIG/config_exp_Sp.ini

## Merge FITS headers
shapepipe_run -c $SP_CONFIG/config_exp_Mh.ini

## Mask images
shapepipe_run -c $SP_CONFIG/config_exp_Ma.ini

## Detect objects
shapepipe_run -c $SP_CONFIG/config_exp_Sx.ini

## Select star candidates
shapepipe_run -c $SP_CONFIG/config_exp_Se.ini

## Validation: Create stats plots
stats_global -o stats -v -c $SP_CONFIG/config_stats.ini

## Create PSF model
shapepipe_run -c $SP_CONFIG/config_exp_Psm.ini

## Interpolate PSF model to star positions (for validation)
shapepipe_run -c $SP_CONFIG/config_exp_Psi.ini

## Validation: PSF residuals
shapepipe_run -c $SP_CONFIG/config_exp_Cp.ini

mkdir psf_validation
shapepipe_run -c $SP_CONFIG/config_exp_Mst.ini
