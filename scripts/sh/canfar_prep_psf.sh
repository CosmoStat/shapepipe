#!/usr/bin/env bash


## Paths
export SP_ROOT=$HOME/astro/repositories/gitlab.cea/ShapePipe

psfval_file_base="validation_psf"
export dir_individual="psf_validation_ind"
export dir_merged="psf_validation_merged"
export fname_merged="psf_cat_full.fits"
pwd=`pwd`


## Functions
function link_s () {
    dst=$1
    src=$2
 
    if [ -e "$src" ]; then
        #echo "link $src already exists, skipping..."
    else
        #echo "create link $dst <- $src"
        ln -s $dst $src
    fi
}


### Start ###

# Create output dirs
for dir in $dir_individual $dir_merged; do
    if [ ! -d "$dir" ]; then
        mkdir -p $dir
    fi
done

# Go through all ngmix output files and extract psf validation files
NGMIX=ngmix*.tgz
for ng in $NGMIX; do
    ID=`echo $ng | perl -ane 's/ngmix_//; s/\.tgz//; print'`
    echo $ID

    valpsf=psfexinterp_$ID.tgz
    tar xf $valpsf
    FILES=output/*/psfexinterp_runner/output/${psfval_file_base}*
    for val in $FILES; do
        base=`basename $val`
        link_s "$pwd/$val" "$dir_individual/$base"
    done
done

# Create merged PSF validation catalog
#$SP_ROOT/scripts/python/merge_star_cat.py -i $dir_individual -o $dir_merged/$fname_merged -v

# Create plots
#$SP_ROOT/scripts/python/MeanShapes.py -o $dir_merged -i $dir_merged/$fname_merged -v -x 20 --max_e=0.05 --max_d=0.005

#tar czf p.tgz psf_validation_merged/*.png
