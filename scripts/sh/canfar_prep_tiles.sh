#!/usr/bin/env bash

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
run_dir="run_combined"
log_path="$pwd/$out_base/log_run_sp_tile.txt"
INPUT="$pwd/$out_base/run_sp_*"
OUTPUT="$pwd/$out_base/$run_dir"
mkdir -p $OUTPUT

# Directories and file patterns to create/link
DIRS=(
	"sextractor_runner"
	"spread_model_runner"
	"psfexinterp_runner"
	"ngmix_runner"
)

PATTERNS=(
	"sexcat_sexcat-*"
	"sexcat_sm-*"
	"galaxy_psf-*"
        "ngmix-*"
)

# Create links
for n in "${!PATTERNS[@]}"; do
    pattern=${PATTERNS[$n]}
    dir=$OUTPUT/${DIRS[$n]}/output
    echo $n $pattern $dir
    mkdir -p $dir
    FILES=(`find $INPUT -name "$pattern"`)
    for file in ${FILES[@]}; do
	target=$file
	link_name=$dir/`basename $file`
	link_s $target $link_name
    done
done

# Update log file
modules=`echo ${DIRS[@]} | tr ' ' ,`
echo "./$out_base/$run_dir $modules" >> $log_path
