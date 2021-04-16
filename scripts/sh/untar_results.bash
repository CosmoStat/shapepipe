#!/usr/bin/env bash

# Name: untar_results.bash
# Description: Untar .tgz files = results of ShapePipe runs
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


NAMES=(
        "final_cat"
        "setools_mask"
        "setools_stat"
        "setools_plot"
        "pipeline_flag"
     )

if [ "$psf" == "psfex" ]; then
  NAMES+=(
          "psfex_interp_exp"
         )
elif [ "$psf" == "mccd" ]; then
  NAMES+=(
          "mccd_fit_val_runner"
         )
fi

# Check number of files
for out in ${NAMES[@]}; do
    echo "$out"
    FILES=${out}_*.tgz
    n_files=${#FILES[@]}
    n_ok=0
    n_fail=0
    for file in $FILES; do
	    tar xf $file
        res=$?
        if [ $res == 0 ]; then
            ((n_ok=n_ok+1))
        else
            ((n_fail=n_fail+1))
        fi
    done
    echo " $n_files tgz files, $n_ok successful untar commands, $n_fail failures"
done
