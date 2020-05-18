#!/usr/bin/env bash

## Paths
remote="vos:cfis/cosmostat/kilbinger/results"
local="."

output=(ngmix psfex psfexinterp setools_plot setools_stat logs)

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
for out in ${output[@]}; do
    cmd="$VCP $remote/$out*.tgz $local"
    echo $cmd
    $cmd
done

# Check number of files
for out in ${output[@]}; do
    n_downl=(`ls -l $local/$out*.tgz | wc`)
    echo "$n_downl '$out' result files downloaded"
done
