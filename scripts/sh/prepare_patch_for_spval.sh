#!/usr/bin/env bash

patch=$1

spdir=$HOME/astro/repositories/github/sp_validation

# Galaxy catalogue
mv ~/psfex/final_cat_${patch}.hdf5 .

# Star catalogue
ln -s $HOME/psfex/$patch/output/run_sp_Ms/merge_starcat_runner/output/full_starcat-0000000.fits

# Parameter file
cp $spdir/notebooks/params.py .

# nb/python script
ln -s $spdir/notebooks/validation.py

# Tile number list
ln -s ~/shapepipe/auxdir/CFIS/tiles_202106/tiles_${patch}.txt 

echo "Change patch in params.py"
