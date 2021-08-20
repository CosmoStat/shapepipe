#!/usr/bin/env bash

# Name: prepare_tiles_for_final.bash
# Description: Create shapepipe run directory with
#              links to all `final_cat` fits files
# Author: Martin Kilbinger <martin.kilbinger@cea.fr
# Package: ShapePipe
# Date: 08/2021
# Version: 0.2

type='final'

## Help string
usage="Usage: $(basename "$0") [OPTIONS]
\n\nOptions:\n
   -h\tthis message\n
   -t, --type TYPE\n
   \tInput type, one in 'final', 'mask', default='$type'\n
"

## Parse command line
while [ $# -gt 0 ]; do
  case "$1" in
    -h)
      echo -ne $usage
      exit 0
      ;;
    -t|--type)
      type="$2"
      shift
      ;;
    *)
      echo -ne usage
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


### Start ###

pwd=`pwd`
out_base="output"
log_path="$pwd/$out_base/log_run_sp.txt"

if [[ "$type" == "final" ]]; then
  run_dir="run_sp_combined"
  INPUT="$pwd/$out_base/run_sp_Mc_*"
  DIRS=(
	  "make_catalog_runner"
  )
  PATTERNS=(
	  "final_cat-*"
  )
elif [[ "$type" == "mask" ]]; then
  run_dir="run_sp_combined_mask"
  INPUT="$pwd/$out_base/run_sp_tile_Ma_*"
  DIRS=(
    "mask_runner"
  )
  PATTERNS=(
    "pipeline_flag-*"
  )
else
  echo "Unknown type '$type'"
  exit 2
fi

OUTPUT="$pwd/$out_base/$run_dir"
mkdir -p $OUTPUT



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
