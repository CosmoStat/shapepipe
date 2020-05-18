#!/usr/bin/env bash

## Paths
local="."

output=(ngmix psfex psfexinterp setools_plot setools_stat logs)


# Check number of files
for out in ${output[@]}; do
    echo "$out"
    FILES=$out*.tgz
    for file in $FILES; do
        echo " $file"
	tar xf $file
	rm $file
    done
done
