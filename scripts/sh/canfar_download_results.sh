#!/usr/bin/env bash

# Name: canfar_download_results.sh
# Description: Download ShapePipe results (.tgz files)
#              from canfar with vos
# Author: Martin Kilbinger <martin.kilbinger@cea.fr>
# Date: 05/2020, 01/2021
# Package: shapepipe


# Command line arguments

## Default values
INPUT_VOS="cfis/cosmostat/kilbinger/results"
INPUT_IDs=""

## Help message
usage="Usage: $(basename "$0") [OPTIONS]\n\n
Options:\n
   -h\tThis message\n
   -i, --input_IDs INPUT_IDs\n
\tASCII file with tile IDs to download,\n
\t default: download all available IDs\n
   --input_vos\n
\tinput path on vos, default=$INPUT_VOS"

## Parse command line
while [ $# -gt 0 ]; do
  case "$1" in
    -h)
      echo -ne $usage
      exit 0
      ;;
    -i|--input_IDs)
      INPUT_IDs="$2"
      shift
      ;;
    --input_vos)
      INPUT_VOS="$2"
      shift
      ;;
    *)
      echo "Invalid option or argument \"$1\""
      exit 2
      ;;
  esac
  shift
done


if [ "$INPUT_IDs" != "" ]; then
    IDs=(`cat $INPUT_IDs`)
    echo "Downloading ${#IDs[@]} IDs"
else
    echo "Downloading all remote IDs"
fi

echo "Remote source = \"$INPUT_VOS\""

# Paths
remote="vos:$INPUT_VOS"
local="."

NAMES=(
        "final_cat"
        "logs"
        "psfex"
        "psfexinterp_exp"
        "setools_mask"
        "setools_stat"
        "setools_plot"
        "pipeline_flag"
     )


# VCP options
# VCP without "vflag" to avoid output to stderr
export VERBOSE=1

if [ $VERBOSE == 1 ]; then
   vflag="-v"
else
   vflag=""
fi
export VCP="vcp $vflag"


### Start ###

# Download files
for name in ${NAMES[@]}; do
    if [ ${#IDs[@]} == 0 ]; then
        cmd="$VCP $remote/$name*.tgz $local"
        $cmd
    else
        for ID in ${IDs[@]}; do
            cmd="$VCP $remote/${name}_$ID.tgz $local"
            $cmd
        done
    fi
done


# TODO: Save IDs in text file, check that
# all files are downloaded correctly for
# each ID.

# Check number of files
for name in ${NAMES[@]}; do
    #n_downl=(`ls -l $local/${name}_*.tgz | wc`)
    set -- `ls $local/${name}_*.tgz 2> /dev/null`
    n_downl="$#"
    echo "$n_downl '$name' result files downloaded"
done
