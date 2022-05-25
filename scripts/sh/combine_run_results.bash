#!/usr/bin/env bash

# Name: combine_run_outputs.bash
# Description: Combine the output results of different ShapePipe runs.
#              Create shapepipe run directory with
#              links to all `final_cat` fits files.
#              In addition, the pipieline run log file is updated.
# Author: Martin Kilbinger <martin.kilbinger@cea.fr>
# Date: 06/2020
# Version: 0.2 (05/2022)

# Command line arguments

## Default values
cat='final'  

## Help string
usage="Usage: $(basename "$0") [OPTIONS]
\n\nOptions:\n
   -h\tthis message\n
   -c, --cat TYPE\n
    \tcatalogue type, one in ['final'|'flag'|'image'], default='$cat'\n
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
  echo "cat (option -c) needs to be in ['final'|'flag'|'image']"
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
elif [ "$cat" == "image" ]; then
  run_dir="run_sp_combined_image"
  INPUT="$pwd/$out_base/run_sp_*Gie_*"
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
elif [ "$cat" == "image" ]; then
  DIRS=(
    "get_images_runner_run_1"
    "get_images_runner_run_2"
  )
  PATTERNS=(
    "CFIS_[iw]*.fits*"
    "[iwf]*-*.fitsfz"
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
modules=`echo ${DIRS[@]} | perl -ane 's/_run_\d{1}//g; print' | tr ' ' ,`
echo "./$out_base/$run_dir $modules" >> $log_path
