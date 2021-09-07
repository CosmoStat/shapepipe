#!/usr/bin/env bash

# Name: post_proc_sp.bash
# Description: Post-process downloaded canfar results and
#              creates final merge catalog.
# Author: Martin Kilbinger <martin.kilbinger@cea.fr>
# Date: 06/2020
# Package: shapepipe

# The shapepipe python virtual environment needs to 
# be active to run this script.


# Command line arguments

## Default values
psf='mccd'
random=0

## Help string
usage="Usage: $(basename "$0") [OPTIONS]
\n\nOptions:\n
   -h\tthis message\n
   -p, --psf MODEL\n
    \tPSF model, one in ['psfex'|'mccd'], default='$psf'\n
   -r, --random\n
   \tcreate random catalogue\n
"

## Parse command line
while [ $# -gt 0 ]; do
  case "$1" in
    -h)
      echo -ne $usage
      exit 0
      ;;
    -p|--psf)
      psf="$2"
      shift
      ;;
    -r|--random)
      random=1
      ;;
    *)
      echo -ne usage
      exit 1
      ;;
  esac
  shift
done


## Check options
if [ "$psf" != "psfex" ] && [ "$psf" != "mccd" ]; then
  echo "PSF (option -p) needs to be 'psfex' or 'mccd'"
  exit 2
fi

# Paths
SP_BASE=$HOME/astro/repositories/github/shapepipe
SP_CONFIG=$SP_BASE/example/cfis


# Merge all psfinterp results and compute PSF residuals
psf_residuals -p $psf

# Prepare output directory with links to all 'final_cat' result files
prepare_tiles_for_final

# Merge final output files to single numpy catalog
input_final=output/run_sp_combined/make_catalog_runner/output
merge_final_cat -i $input_final -p $SP_CONFIG/final_cat.param -v 

if [ "$random" == "1" ]; then
  # Create random catalogues for each tile
  shapepipe_run -c $SP_CONFIG/config_Rc.ini

  # Merge all random catalogues to single numpy catalogue
  input_final=output/run_sp_Rc/random_cat_runner/output
  merge_final_cat -i $input_final -n random_cat -v 
fi

# Merge star catalogue and plot PSF residuals
shapepipe_run -c $SP_CONFIG/config_MsPl_mccd.ini 
