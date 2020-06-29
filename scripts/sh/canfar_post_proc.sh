#!/usr/bin/env bash

# Name: canfar_post_proc.sh
# Description: Post-process downloaded canfar results and
#              creates final merge catalog.
# Author: Martin Kilbinger <martin.kilbinger@cea.fr>
# Date: 06/2020
# Package: shapepipe


# The shapepipe python virtual environment needs to 
# be active to run this script.

SP_BASE=$HOME/astro/repositories/github/shapepipe
SP_CONFIG=$SP_BASE/example/cfis

# To download results from canfar, use
#
# canfar_download_results.sh
#
# On candide this needs to be done on
# the login node.

# Un-tar all .tgz results files
$SP_BASE/scripts/sh/untar_results.sh

# Merge all psfinterp results and compute PSF residuals
$SP_BASE/scripts/sh/canfar_psf_residuals.sh

# Merge final output files to single mother catalog
input_final=output/run_sp_combined/make_catalog_runner/output
$SP_BASE/scripts/python/merge_final_cat.py -i $input_final -p $SP_CONFIG/final_cat.param -v 
