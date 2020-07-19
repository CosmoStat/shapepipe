#!/usr/bin/env bash

# Name: canfar_prep_tiles.sh
# Description: Create shapepipe run directory with
#              links to all `final_cat` fits files
# Author: Martin Kilbinger <martin.kilbinger@cea.fr
# Package: ShapePipe
# Date: 06/2020
# Version: 0.1

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


### Start ###

pwd=`pwd`
out_base="output"
run_dir="run_sp_combined"
log_path="$pwd/$out_base/log_run_sp_tile.txt"
INPUT="$pwd/$out_base/run_sp_Mc_*"
OUTPUT="$pwd/$out_base/$run_dir"
mkdir -p $OUTPUT

# Directories and file patterns to create/link
DIRS=(
	"make_catalog_runner"
)

PATTERNS=(
	"final_cat-*"
)

# Create links
for n in "${!PATTERNS[@]}"; do
    pattern=${PATTERNS[$n]}
    dir=$OUTPUT/${DIRS[$n]}/output
    echo $n $pattern $dir
    mkdir -p $dir
    FILES=(`find $INPUT -name "$pattern"`)
    n_files=${#FILES[@]}
    i=0
    for file in ${FILES[@]}; do
	    target=$file
	    link_name=$dir/`basename $file`
	    link_s $target $link_name
        ((i=i+1))
    done
    echo " $n_files target files, $i links created/skipped"
done

# Update log file
modules=`echo ${DIRS[@]} | tr ' ' ,`
echo "./$out_base/$run_dir $modules" >> $log_path
