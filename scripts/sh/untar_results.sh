#!/usr/bin/env bash

# Description: Untar .tgz files = results of ShapePipe runs
# Author: Martin Kilbinger <martin.kilbinger@cea.fr>
# Date: 05/2020
# Package: shapepipe


## Paths
local="."

output=(ngmix psfex psfexinterp_exp psfexinterp_tile setools_plot setools_stat setools_mask logs)
NAMES=(
        "psfex"
        "psfexinterp_exp"
        "setools_mask"
        "setools_stat"
        "setools_plot"
        "sextractor"
        "spread_model"
        "psfexinterp_tile"
	"ngmix"
	"logs"
     )



# Check number of files
for out in ${NAMES[@]}; do
    echo "$out"
    FILES=${out}_*.tgz
    for file in $FILES; do
        echo " $file"
	tar xf $file
	#rm $file
    done
done
