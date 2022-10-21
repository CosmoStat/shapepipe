#!/usr/bin/env bash

# Name: prepare_star_cat.bash
# Description: Create directory and links to all PSF or star catalogue files
#              from previous ShapePipe runs.
# Author: Martin Kilbinger <martin.kilbinger@cea.fr>


# Command line arguments

## Default values
psf='mccd'

## Help string
usage="Usage: $(basename "$0") [OPTIONS]
\n\nOptions:\n
   -h\tthis message\n
   -p, --psf MODEL\n
    \tPSF model, one in ['psfex'|'mccd'|'setools'], default='$psf'\n
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


## Path variables
if [ "$psf" == "psfex" ] || [ "$psf" == "mccd" ]; then
  psfval_file_base="validation_psf"
  dir_individual="psf_validation_ind"
else
  psfval_file_base="mask/star_selection"
  dir_individual="star_all_ind"
fi

pwd=`pwd`


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


# Create output dirs
if [ ! -d "$dir_individual" ]; then
    mkdir -p $dir_individual
fi

if [ "$psf" == "psfex" ]; then
  runner="psfex_interp_runner"
elif [ "$psf" == "mccd" ]; then
  runner="mccd_fit_val_runner"
else
  runner="setools_runner"
fi

# Find all psf validation files and create links.
# Assumes untar_results.sh has been run before.
n_skipped=0
n_created=0
FILES=output/*/${runner}/output/${psfval_file_base}*
for val in ${FILES[@]}; do
    base=`basename $val`
    link_s "$pwd/$val" "$dir_individual/$base"
done
echo " Created $n_created links, skipped $n_skipped files"
