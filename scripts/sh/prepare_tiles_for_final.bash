#!/usr/bin/env bash

# Name: prepare_tiles_for_final.bash
# Description: Create shapepipe run directory with
#              links to all `final_cat` fits files
# Author: Martin Kilbinger <martin.kilbinger@cea.fr>
# Date: 06/2020
# Version: 0.1

# Command line arguments

## Default values
cat='final'  

## Help string
usage="Usage: $(basename "$0") [OPTIONS]
\n\nOptions:\n
   -h\tthis message\n
   -c, --cat TYPE\n
    \tCatalogue type, one in ['final'|'flag'|'image'], default='$cat'\n
"

## Parse command line
while [ $# -gt 0 ]; do
  case "$1" in
    -h)
      echo -ne $usage
      exit 0
      ;;
    -c|--cat)
      cat="$2"
      shift
      ;;
    *)
      echo -ne $usage
      exit 1
      ;;
  esac
  shift
done


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

## Check options
if [ "$cat" != "final" ] && [ "$cat" != "flag" ] && [ "$cat" != "image" ]; then
  echo "cat (option -c) needs to be 'final', 'flag', or 'image'"
  exit 2
fi


### Start ###

pwd=`pwd`
out_base="output"

if [ "$cat" == "final" ]; then
  run_dir="run_sp_combined"
  INPUT="$pwd/$out_base/run_sp_Mc_*"
elif [ "$cat" == "flag" ]; then
  run_dir="run_sp_combined_flag"
  INPUT="$pwd/$out_base/run_sp_tile_Ma_*"
else
  run_dir="run_sp_combined_image"
  INPUT="$pwd/$out_base/run_sp_Git_*"
fi

log_path="$pwd/$out_base/log_run_sp.txt"
OUTPUT="$pwd/$out_base/$run_dir"
mkdir -p $OUTPUT

# Directories and file patterns to create/link
if [ "$cat" == "final" ]; then
  DIRS=(
	  "make_catalog_runner"
  )
  PATTERNS=(
	  "final_cat-*"
  )
elif [ "$cat" == "flag" ]; then
  DIRS=(
	  "mask_runner"
  )
  PATTERNS=(
	  "pipeline_flag-*"
  )
else
  DIRS=(
    "get_images_runner"
  )
  PATTERNS=(
	  "CFIS_image-*"
  )
fi

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
