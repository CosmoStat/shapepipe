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

# Merge seperated ngmix catalogs
shapepipe_run -c $SP_CONFIG/config_merge_sep_cats.ini 

# Create links to required tiles results
$SP_BASE/scripts/sh/canfar_prep_tiles.sh

# Create catalogs with all necessary tile information
shapepipe_run -c $SP_CONFIG/config_make_cat.ini

# Merge final output files to single mother catalog
$SP_BASE/scripts/python/merge_final_cat.py
