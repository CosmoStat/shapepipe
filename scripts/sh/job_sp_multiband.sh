#!/usr/bin/env bash

# Name: job_sp_multiband.sh
# Description: Process UNIONS tiles to
#              create multi-band catalogue
# Author: Martin Kilbinger, <martin.kilbinger@cea.fr>
# Date: 13/11/2020
# Package: ShapePipe


if [ $# == 0 ]; then
    echo "Usage:"
    echo "  `basename $0` TILE_ID_1 [TILE_ID_2 [...]]"
    echo "    TILE_ID = xxx.yyy"
    echo "Examples:"
    echo "  `basename $0` 282.247"
    echo "  `basename $0` 282.247 282.238"
    exit 1
fi

# Set-up

# Path variables
export SP_HOME=$HOME/astro/repositories/github/shapepipe
export SP_CONFIG=$SP_HOME/example/unions
export SP_RUN=.

export LD_LIBRARY_PATH=$HOME/.conda/envs/shapepipe/lib

# Command line arguments
TILE_ARR=($@)
n_tile=${#TILE_ARR[@]}
echo "Processing $n_tile tile(s)"

# Create input and output directories
mkdir -p $SP_RUN
cd $SP_RUN
mkdir -p output

# Write tile numbers to ASCII input file
rm -rf tile_numbers.txt
for TILE in ${TILE_ARR[@]}; do
   echo $TILE >> tile_numbers.txt
done

# Processing: run pipeline

# Download
shapepipe_run -c $SP_CONFIG/config_get_tiles_symlink.ini
#shapepipe_run -c $SP_CONFIG/config_get_vcp.ini

# Convert Pan-STARRS image names
#$SP_HOME/scripts/python/ps_convert_file_names.py

# Unzip FITS files and/or remove unused HDU:
# CFIS weights; UNIONS images, weights, flags
shapepipe_run -c $SP_CONFIG/config_unfz.ini

# Create flags for CFIS r-band images: add star, halo, and Messier
# object masks.
shapepipe_run -c $SP_CONFIG/config_tile_mask_r.ini

# Detect objects on r-band images, measure properties
# on u-, r-, i-, z-band images
shapepipe_run -c $SP_CONFIG/config_tile_detect_r.ini
shapepipe_run -c $SP_CONFIG/config_tile_detect_u.ini
shapepipe_run -c $SP_CONFIG/config_tile_detect_i.ini
shapepipe_run -c $SP_CONFIG/config_tile_detect_z.ini

# Merge catalogs
shapepipe_run -c $SP_CONFIG/config_tile_paste_cat.ini
