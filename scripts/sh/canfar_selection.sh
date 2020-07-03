#!/usr/bin/env bash

# Name: canfar_selection.sh
# Description: Create selection of canfar-run tiles for
#              post-processing.
# Author: Martin Kilbinger <martin.kilbinger@cea.fr>
# Date: 06/2020
# Package: shapepipe


# Command line arguments
if [ $# != 2 ] && [ $# != 3 ]; then
    echo "Usage: canfar_selection.sh IN_DIR LIST_PATH [OUT_DIR]"
    echo "  IN_DIR: input directory with .tgz files"
    echo "  LIST_PATH: ASCII file with tile IDs"
    echo "  OUT_DIR: output directory name, default=tile ID file base name"
    exit 1
fi

IN_DIR=$1
LIST_PATH=$2

if [ $# ==  3 ]; then
    OUT_DIR=$3
else
    OUT_DIR=`echo "${LIST_PATH%.*}"`
fi

echo "IN_DIR=$IN_DIR"
echo "LIST_PATH=$LIST_PATH"
echo "OUT_DIR=$OUT_DIR"


## Functions
function link_s () {
    target=$1
    link_name=$2

    if [ -e "$link_name" ]; then
        echo "link with name $link_name already exists, skipping..."
    else
        echo "create link $target <- $link_name"
        ln -s $target $link_name
    fi
}


# Variables

pwd=`pwd`

# Tar archive base names
NAMES=(
        "psfex"
        "psfexinterp_exp"
        "setools_mask"
        "setools_stat"
        "setools_plot"
        "final_cat"
	    "logs"
     )

## Start
if [[ -d $OUT_DIR ]]; then
    echo "Error: Output directory $OUT_DIR already exists, exiting..."
    exit 2
fi

mkdir $OUT_DIR
for i in `cat $LIST_PATH`; do
    for out in ${NAMES[@]}; do
        target=$pwd/$IN_DIR/${out}_$i.tgz
        link_name=$OUT_DIR/${out}_$i.tgz
        link_s $target $link_name
    done
done
