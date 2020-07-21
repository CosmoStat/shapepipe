#!/usr/bin/env bash

# Name: canfar_download_results.sh
# Description: Download ShapePipe results (.tgz files)
#              from canfar with vos
# Author: Martin Kilbinger <martin.kilbinger@cea.fr>
# Date: 05/2020
# Package: shapepipe

if [ "$1" == "-h" ]; then
    echo "Usage:"
    echo "  canfar_download_all.sh [-h] [ID_FILE]"
    echo "     ID_FILE     ASCII file with tile IDs to download,"
    echo "                  default: download all available IDs"
    echo "     -h          This message"
    exit 0
else
    if [ "$#" == "1" ]; then
        IDs=(`cat $1`)
        echo "Downloading ${#IDs[@]} IDs"
    fi
fi

## Paths
#remote="vos:cfis/cosmostat/kilbinger/results_lsb"
remote="vos:cfis/cosmostat/kilbinger/results"
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
export VCP_QUICK=0
export VERBOSE=1

if [ $VCP_QUICK == 1 ]; then
   qflag="--quick"
else
   qflag=""
fi
if [ $VERBOSE == 1 ]; then
   vflag="-v"
else
   vflag=""
fi
export VCP="vcp $qflag $vflag"


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
    echo "$n_downl '$name' result files downloaded"
done
