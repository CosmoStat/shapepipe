#!/usr/bin/env bash

# Name: untar_results.ch
# Description: Untar .tgz files = results of ShapePipe runs
# Author: Martin Kilbinger <martin.kilbinger@cea.fr>
# Date: 05/2020
# Package: shapepipe


NAMES=(
        "final_cat"
        "psfex"
        "psfex_interp_exp"
        "setools_mask"
        "setools_stat"
        "setools_plot"
        "pipeline_flag"
     )


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
