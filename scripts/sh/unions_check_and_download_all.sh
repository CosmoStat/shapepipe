#!/usr/bin/env bash

# Name: unions_check_and_download_all.bash
# Description: Compare number of local and remote UNIONS images
#              and download missing files.
#              Run this script on candide.
# Author: Martin Kilbinger <martin.kilbinger@cea.fr>
# Date: 01/2021
# Package: ShapePipe

# Global variables

remote_base="vos:cfis"

# Functions

## Check numbers of remote and local files
function numbers () {
   remote=$1
   dir=$2
   pattern=$3

   for p in ${pattern[@]}; do
      printf "%35s %20s %50s " $remote $dir $p

      cmd="ls -l $dir/$p"
      nloc=(`$cmd | wc`)
      printf "%8s " ${nloc[0]}

      cmd="vls $remote/$p" 
      nrem=(`$cmd | wc`)
      printf "%8s" ${nrem[0]}

      printf "\n"
   done
}

# Start

cd ~/astro/data/CFIS

if [[ "$CONDA_DEFAULT_ENV" != "shapepipe" ]]; then
   conda activate shapepipe
fi

# Output header
printf "%35s %20s %50s %8s %8s\n" "remote" "dir" "pattern" "nloc" "nrem"

dir="tiles_DR2"
remote="$remote_base/$dir"
pattern=("CFIS.???.???.r.fits" "CFIS.???.???.r.weight.fits.fz" "CFIS.???.???.u.fits" "CFIS.???.???.u.weight.fits.fz")
#pattern=("CFIS.005.253.r.fits" "CFIS.005.253.r.weight.fits.fz" "CFIS.005.253.u.fits" "CFIS.005.253.u.weight.fits.fz")
numbers $remote $dir $pattern

dir="Pan-STARRS/ps_orig"
remote="$remote_base/panstarrs/DR2/skycell.???"
#remote="$remote_base/panstarrs/DR2/skycell.005"
pattern=("CFIS.V0.skycell.???.???.stk.*.unconv.fits" "CFIS.V0.skycell.???.???.stk.*.unconv.wt.fits" "CFIS.V0.skycell.???.???.stk.*.unconv.mask.fits")
#pattern=("CFIS.V0.skycell.005.253.stk.*.unconv.fits" "CFIS.V0.skycell.005.253.stk.*.unconv.wt.fits" "CFIS.V0.skycell.005.253.stk.*.unconv.mask.fits")
numbers $remote $dir $pattern

basedir="pitcairn"
dir="${basedir}_DR2"
remote="$remote_base/$basedir"
pattern=("2??????p.fits.fz")
numbers $remote $dir $pattern

basedir="weights"
dir="${basedir}_DR2"
remote="$remote_base/$basedir"
pattern=("2??????p.weight.fits.fz")
numbers $remote $dir $pattern

basedir="flags"
dir="${basedir}_DR2"
remote="$remote_base/$basedir"
pattern=("2??????p.flag.fits.fz")
numbers $remote $dir $pattern
