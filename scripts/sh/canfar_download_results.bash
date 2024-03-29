#!/usr/bin/env bash

# Name: canfar_download_results.bash
# Description: Download ShapePipe results (.tgz files)
#              from canfar with vos
# Author: Martin Kilbinger <martin.kilbinger@cea.fr>
# Date: v1.0 05/2020
#       v1.1 01/2021

# Command line

## Default parameters
INPUT_VOS="cosmostat/kilbinger/results"
VERBOSE=0
psf="mccd"
only_mask=0


usage="Usage: $(basename "$0") [OPTIONS]
\n\nOptions:\n
    -h\tthis message\n
    -i, --input_IDs ID_FILE\n
    \tASCII file with tile IDs to download, default:\n
    \tdownload all available IDs\n
    --input_vos PATH\n
    \tinput path on vos:cfis, default='$INPUT_VOS'\n
  -p, --psf MODEL\n
    \tPSF model, one in ['psfex'|'mccd'], default='$psf'\n
  -m\tonly mask files\n
  -v\tverbose output\n
"
  
## Parse command line
while [ $# -gt 0 ]; do
  case "$1" in
    -h)
      echo -ne $usage
      exit 0
      ;;
    -i|--input_IDs)
      IDs=(`cat $2`)
      echo "Downloading ${#IDs[@]} ID(s)"
      shift
      ;;
    --input_vos)
      INPUT_VOS="$2"
      shift
      ;;
    -p|--psf)
      psf="$2"
      shift
      ;;
    -m)
      only_mask=1
      ;;
    -v)
      VERBOSE=1
      ;;
    *)
      echo "Invalid command line argument '$1'"
      echo -ne $usage
      exit 1
      ;;
  esac
  shift
done

## Check options
if [ "$psf" != "psfex" ] && [ "$psf" != "mccd" ]; then
  echo "PSF (option -p) needs to be 'psf' or 'mccd'"
  exit 2
fi

## Paths
remote="vos:cfis/$INPUT_VOS"
local="."

if [ $only_mask == 1 ]; then
  NAMES=("pipeline_flag")
else
  NAMES=(
        "final_cat"
        "logs"
        "setools_mask"
        "setools_stat"
        "setools_plot"
        "pipeline_flag"
  )

  if [ $psf == "psfex" ]; then
    NAMES+=(
          "psfex_interp_exp"
    )
  else
    NAMES+=(
          "mccd_fit_val_runner"
    )
  fi
fi

if [ $VERBOSE == 1 ]; then
   vflag="-v"
else
   vflag=""
fi
export VCP="vcp $vflag"


### Start ###

# Download files
for name in ${NAMES[@]}; do
    if [ ${#IDs[@]} == 0 ]; then
        cmd="$VCP $remote/$name*.tgz $local"
        $cmd
    else
        for ID in ${IDs[@]}; do
            cmd="$VCP $remote/${name}_$ID.tgz $local"
            $cmd
        done
    fi
done

# Check number of files
for name in ${NAMES[@]}; do
    n_downl=(`ls -l $local/${name}_*.tgz | wc`)
    echo "$n_downl '$name' result files downloaded from $remote"
done
