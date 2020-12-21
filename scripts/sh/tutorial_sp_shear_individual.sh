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


# Functions

## Print string, executes command, and prints return value.
function command () {
   cmd=$1
   str=$2

   echo "$str: running '$cmd'"
   $cmd
   res=$?
   if [ $res != 0 ]; then
     echo -e "error, return value = $res"
     if [ $STOP == 1 ]; then
       echo "exiting $(basename "$0"), error in command '$cmd'"
       exit $res
     else
       echo "continuing $(basename "$0"), error in command '$cmd'"
     fi
   fi

   return $res
}

## Run shapepipe command
command_sp() {
   cmd=$1
   str=$2

   STOP=0
   command "$1" "$2"
   res=$?
   if [ $res != 0 ]; then
      echo "exiting $(basename "$0"), '$cmd' returned $res"
      exit $res
   fi

}


# Command line arguments

## Parse command line
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

# Work-around of SExtractor&psfex lib bug
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib

# Create output path
mkdir -p output

# Write tile numbers to ASCII input file
rm -rf tile_numbers.txt
for id in ${ID[@]}; do
   echo $id >> tile_numbers.txt
done


# Pre-processing

## Select tile IDs to process
#cfis_field_select -i $SP_BASE/aux/CFIS/tiles_202007/tiles_all_order.txt --coord 213.68deg_54.79deg -t tile --input_format ID_only --out_name_only --out_ID_only -s -o tile_numbers -v

## Retrieve tiles
if [ "$retrieve" == "vos" ]; then

  command_sp "shapepipe_run -c $SP_CONFIG/config_tile_Gi.ini" "Get tiles ($retrieve)"

elif [ "$retrieve" == "symlink" ]; then

  if [ ! -d "data" ]; then
    echo "Input directory 'data' not found"
    exit 2
  fi
  command_sp "shapepipe_run -c $SP_CONFIG/config_tile_Gi_symlink.ini" "Get tiles ($retrieve)"

else

   echo "Invalid retrieve value '$retrieve'"

fi

command_sp "shapepipe_run -c $SP_CONFIG/config_tile_Uz.ini" "Uncompress tile weights"

command_sp "shapepipe_run -c $SP_CONFIG/config_tile_Fe.ini" "Find exposures"

## Retrieve single exposures
if [ "$retrieve" == "vos" ]; then

  command_sp "shapepipe_run -c $SP_CONFIG/config_exp_Gi.ini" "Get exposures ($retrieve)"

elif [ "$retrieve" == "symlink" ]; then

  command_sp "shapepipe_run -c $SP_CONFIG/config_exp_Gi_symlink.ini" "Get exposures ($retrieve)"

fi

# Processing of single exposures

command_sp "shapepipe_run -c $SP_CONFIG/config_exp_Sp.ini" "Exp: split into single-HDU files"

command_sp "shapepipe_run -c $SP_CONFIG/config_exp_Mh.ini" "Exp: merge headersss"

command_sp "shapepipe_run -c $SP_CONFIG/config_exp_Ma.ini" "Exp: mask"

command_sp "shapepipe_run -c $SP_CONFIG/config_exp_Sx.ini" "Exp: detect star candidates"

command_sp "shapepipe_run -c $SP_CONFIG/config_exp_Se.ini" "Exp: select stars"

stats_global -o stats -v -c $SP_CONFIG/config_stats.ini

command_sp "shapepipe_run -c $SP_CONFIG/config_exp_Psm.ini" "Exp: Create PSF model"

command_sp "shapepipe_run -c $SP_CONFIG/config_exp_Psi.ini" "Exp: Interpolate PSF to star positions"

## Validation: PSF residuals

command_sp "shapepipe_run -c $SP_CONFIG/config_exp_Cp.ini" "Exp: Combing PSF runs"

### Merge PSF files
mkdir psf_validation
command_sp "shapepipe_run -c $SP_CONFIG/config_exp_Mst.ini" "Exp: Merge PSF files"

### Create plots
MeanShapes -o psf_validation -i psf_validation/full_starcat.fits -v -x 20

# Some bad hacks to get additional input files...
input_psfex=`find . -name star_split_ratio_80-*.psf | head -n 1`
ln -s `dirname $input_psfex` input_psfex


# Processing of tiles

command_sp "shapepipe_run -c $SP_CONFIG/config_tile_Ma.ini" "Tile: mask"

command_sp "shapepipe_run -c $SP_CONFIG/config_tile_Sx.ini" "Tile: Detect galaxy candidates"
