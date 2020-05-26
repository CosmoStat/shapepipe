#!/usr/bin/env bash

## Paths
remote="vos:cfis/cosmostat/kilbinger/results"
local="."

output=(ngmix psfex psfexinterp setools_plot setools_stat logs)
NAMES=(
        "psfex"
        "psfexinterp_exp"
        "setools_mask"
        "setools_stat"
        "setools_plot"
        "sextractor"
        "spread_model"
        "psfexinterp_tile"
     )
NAMES=( "${NAMES[@]}" "ngmix" "logs")


# VCP options
# VCP without "vflag" to avoid output to stderr
export VCP_QUICK=1
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
    cmd="$VCP $remote/$name*.tgz $local"
    echo $cmd
    $cmd
done

# TODO: Save IDs in text file, check that
# all files are downloaded correctly for
# each ID.

# Check number of files
for name in ${NAMES[@]}; do
    n_downl=(`ls -l $local/$name_*.tgz | wc`)
    echo "$n_downl '$name' result files downloaded"
done
