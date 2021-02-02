#!/usr/bin/env bash

# Description: Compute and plot PSF residuals from
#	       results processed on canfar
# Author: Martin Kilbinger <martin.kilbinger@cea.fr>
# Date: 05/2020
# Package: shapepipe

# Command line arguments

## Default values
psf='mccd'

## Help string
usage="Usage: $(basename "$0") [OPTIONS]
\n\nOptions:\n
   -h\tthis message\n
   -p, --psf MODEL\n
    \tPSF model, one in ['psfex'|'mccd'], default='$psf'\n
"

## Parse command line
while [ $# -gt 0 ]; do
  case "$1" in
    -h)
      echo -ne $usage
      exit 0
      ;;
    -p|--psf)
      psf="$2"
      shift
      ;;
    *)
      echo -ne usage
      exit 1
      ;;
  esac
  shift
done


## Default variables
psfval_file_base="validation_psf"
dir_individual="psf_validation_ind"
dir_merged="psf_validation_merged"
fname_merged="psf_cat_full.fits"
pwd=`pwd`

# Command line sarguments

## Help string
usage="Usage: $(basename "$0") [OPTIONS]
\n\nOptions:\n
   -h\tthis message\n
   -p, --psf MODEL\n
    \tPSF model, one in ['psfex'|'mccd'], default='$psf'\n
"

## Parse command line
TILE_ARR=()
while [ $# -gt 0 ]; do
  case "$1" in
    -h)
      echo -ne $usage
      exit 0
      ;;
    -p|--psf)
      psf="$2"
      shift
      ;;
    *)
      echo "Invalid command line argument '$2'"
      echo -ne $usage
      exit 1
      ;;
  esac
  shift
done

## Functions
function link_s () {
    target=$1
    link_name=$2
 
    if [ -L "$link_name" ]; then
	      let "n_skipped+=1"
    else
        ln -s $target $link_name
	      let "n_created+=1"
    fi

    return $n
}


### Start ###

# Create output dirs
for dir in $dir_individual $dir_merged; do
    if [ ! -d "$dir" ]; then
        mkdir -p $dir
    fi
done

if [ "$psf" == "psfex" ]; then
  runner="psf_interp_runner"
else
  runner="mccd_fit_val_runner"
fi

# Find all psf validation files and create links.
# Assumes untar_results.sh has been run before.
n_skipped=0
n_created=0
FILES=output/*/${runner}/output/${psfval_file_base}*
for val in $FILES; do
    base=`basename $val`
    link_s "$pwd/$val" "$dir_individual/$base"
done
echo " Created $n_created links, skipped $n_skipped files"

# Create merged PSF validation catalog
merge_star_cat_${psf} -i $dir_individual -o $dir_merged/$fname_merged -v

# Create plots
MeanShapes -o $dir_merged -i $dir_merged/$fname_merged -v -x 20 --max_e=0.05 --max_d=0.005

#tar czf p.tgz psf_validation_merged/*.png

### End ###
