#!/usr/bin/env bash

# Name: canfar_download_results.sh
# Description: Download ShapePipe results (.tgz files)
#              from canfar with vos
# Author: Martin Kilbinger <martin.kilbinger@cea.fr>
# Date: v1.0 05/2020
#       v1.1 01/2021
# Package: shapepipe

# Command line

## Default parameters
INPUT_VOS="cosmostat/kilbinger/results"
VERBOSE=0

usage="Usage: $(basename "$0") [OPTIONS]
\n\nOptions:\n
    -h\tthis message\n
    -i, --input_IDs ID_FILE\n
    \tASCII file with tile IDs to download, default:\n
    \t download all available IDs\n
    --input_vos PATH\n
    \tinput path on vos:cfis, default='$INPUT_VOS'\n
    -v\tverbose output\n
"
  
## Parse command line
while [ $# -gt 0 ]; do
  case "$1" in
    -h)
      echo -ne $usage
      exit 0
      ;;
    -i|--input_IDs)
      IDs=(`cat $2`)
      echo "Downloading ${#IDs[@]} ID(s)"
      shift
      ;;
    --input_vos)
      INPUT_VOS="$2"
      shift
      ;;
    -v)
      VERBOSE=1
      ;;
    *)
      echo "Invalid command line argument '$1'"
      echo -ne $usage
      exit 1
      ;;
  esac
  shift
done


## Paths
remote="vos:cfis/$INPUT_VOS"
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
    n_downl=(`ls -l $local/$name_*.tgz | wc`)
    echo "$n_downl '$name' result files downloaded from $RESULTS"
done
