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

## Help string
usage="Usage: $(basename "$0") [OPTIONS]
\n\nOptions:\n
   -h\tthis message\n
   -p, --psf MODEL\n
    \tPSF model, one in ['psfex'|'mccd'], default='$psf'\n
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
    *)
      echo -ne usage
      exit 1
      ;;
  esac
  shift
done


# Paths
SP_BASE=$HOME/astro/repositories/github/shapepipe
SP_CONFIG=$SP_BASE/example/cfis


# To download results from canfar, use
#
# canfar_download_results.sh
#
# On candide this needs to be done on
# the login node.

# To Un-tar all .tgz results files, use
#
# $SP_BASE/scripts/sh/untar_results.sh
#

# Merge all psfinterp results and compute PSF residuals
psf_residuals -p $psf

# Prepare output directory with links to all 'final_cat' result files
prepare_tiles_for_final.sh

# Merge final output files to single mother catalog
input_final=output/run_sp_combined/make_catalog_runner/output
merge_final_cat -i $input_final -p $SP_CONFIG/final_cat.param -v 
