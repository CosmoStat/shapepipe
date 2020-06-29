#!/usr/bin/env bash

# Name: untar_results.ch
# Description: Untar .tgz files = results of ShapePipe runs
# Author: Martin Kilbinger <martin.kilbinger@cea.fr>
# Date: 05/2020
# Package: shapepipe


NAMES=(
        "psfex"
        "psfexinterp_exp"
        "setools_mask"
        "setools_stat"
        "setools_plot"
        "final_cat"
	    "logs"
     )


# Check number of files
for out in ${NAMES[@]}; do
    echo "$out"
    FILES=${out}_*.tgz
    for file in $FILES; do
        echo " $file"
	tar xf $file
    done
done
