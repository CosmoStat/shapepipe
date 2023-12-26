#!/usr/bin/env bash

# Name: combine_runs.bash
# Description: Create new shapepipe run directory with
#              links to source files from combined existing runs
# Author: Martin Kilbinger <martin.kilbinger@cea.fr>


# Command line arguments

## Default values
cat='final'
psf="mccd"

## Help string
usage="Usage: $(basename "$0") [OPTIONS]
\n\nOptions:\n
   -h\tthis message\n
   -p, --psf MODEL\n
    \tPSF model, one in ['psfex'|'mccd'|'setools'], default='$psf'\n
   -c, --cat TYPE\n
    \tCatalogue type, one in ['final'|'flag'|'psf'|'image'], default='$cat'\n
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
    -p|--psf)
      psf="$2"
      shift
      ;;
    *)
      echo -ne $usage
      exit 1
      ;;
  esac
  shift
done


## Check options
if [ "$cat" != "final" ] \
  && [ "$cat" != "flag" ] \
  && [ "$cat" != "psf" ] \
  && [ "$cat" != "image" ]; then
  echo "cat (option -c) needs to be 'final', 'flag', 'psf', or 'image'"
  exit 2
fi

## Check options
if [ "$psf" != "psfex" ] \
  && [ "$psf" != "mccd" ] \
  && [ "$psf" != "setools" ]; then
  echo "PSF (option -p) needs to be 'psfex' or 'mccd'"
  exit 2
fi


## Functions
function link_s () {
    target=$1
    link_name=$2

    if [ -L "$link_name" ]; then
        echo "link with name $link_name already exists, skipping..."
        let "n_skipped+=1"
    else
        echo "create link $target <- $link_name"
        ln -s $target $link_name
        let "n_created+=1"
    fi
}


# Start program

n_skipped=0
n_created=0

pwd=`pwd`
out_base="output"

# Set paths:
## run_out: target output new run directory
## run_in: source input run base directory
## module: source input module runner sub-directory
## pattern: source file pattern

run_out="run_sp_combined_$cat"

if [ "$cat" == "final" ]; then

  # v1
  #run_in="$pwd/$out_base/run_sp_Mc_*"
  # v2
  run_in="$pwd/tile_runs/*/$out_base/run_sp_Mc_*"

  module="make_catalog_runner"
  pattern="final_cat-*"

elif [ "$cat" == "flag" ]; then

  run_in="$pwd/$out_base/run_sp_MaMa_*/mask_runner_run_1"
  module="mask_runner_run_1"
  pattenr="pipeline_flag-*"

elif [ "$cat" == "image" ]; then

  run_in="$pwd/$out_base/run_sp_Git_*"
  module="get_images_runner"
  pattern="CFIS_image-*"

elif [ "$cat" == "psf" ]; then

  #MKDEBUG TODO: add option
  # v1
  #run_in="$pwf/$out_base/run_sp_exp_Pi_*"
  # v2
  run_in="$pwd/exp_runs/*/$out_base/run_sp_exp_Pi_*"

  pattern="validation_psf-*"
  if [ "$psf" == "psfex" ]; then
    module="psfex_interp_runner"
  elif [ "$psf" == "setools" ]; then
    module="setools_runner"
  else
    module="mccd_interp_runner"
  fi

else

  echo "Invalid catalogue type $cat"
  exit 2

fi


OUTPUT="$pwd/$out_base/$run_out"
mkdir -p $OUTPUT


# Create links

## target directory
outdir=$OUTPUT/$module/output
mkdir -p $outdir

## identify source files

# The following can result in an "Argument list too long" error
#FILES=(`find $run_in -type f -name "$pattern" -print0 | xargs -0 echo`)

i=0
for dir in $run_in; do
  FILES=(`find $dir -type f -name "$pattern" -print0 | xargs -0 echo`)

  ## Look over source files
  for file in ${FILES[@]}; do

  target=$file
  link_name=$outdir/`basename $file`
  link_s $target $link_name
  ((i=i+1))

  done

done

#echo " $n_files target files, $i links created/skipped"
echo " $i total, "n_skipped skipped, "n_created links created"

# Update log file
update_runs_log_file.py
