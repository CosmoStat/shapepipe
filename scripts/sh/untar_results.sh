#!/usr/bin/env bash

## Paths
local="."

output=(ngmix psfex psfexinterp_exp psfexinterp_tile setools_plot setools_stat setools_mask logs)


# Check number of files
for out in ${output[@]}; do
    echo "$out"
    FILES=${out}_*.tgz
    for file in $FILES; do
        echo " $file"
	tar xf $file
	#rm $file
    done
done
