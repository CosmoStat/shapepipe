#!/usr/bin/env bash

# Name: find_exposures.sh
# Description: Script to convert THELI images to shapepipe 
#              bookkeeping conventions
# Author: Lucie Baumont <lucie.baumont@cea.fr>
# Date: v1.0 1/2023


# Command line arguments

## Help string
usage="Usage: $(basename "$0") THELI_DIR OUTPUT_DIR\n
   THELI_DIR\n
   \thome director of THELI images, eg. /n17data/baumont/THELI_W3/W3/CFISCOLLAB_V2.0.0A\n
   OUTPUT_DIR\n
   \troot directory where outputs will go \n
"

## Help if no arguments
if [ -z $1 ]; then
        echo -ne $usage
        exit 1
fi
theli_dir=$1
output_dir=$2

#make necessary directories
mkdir -p ${output_dir}/tiles
mkdir -p ${output_dir}/tiles/exposure_lists
mkdir -p ${output_dir}/exposures
mkdir -p ${output_dir}/exposures/headers

for dir in ${theli_dir}/CFIS*;
    do
    # create link for all tiles and rename them
    coadd_dir=$dir/r.MP9602/coadd_V2.0.0A/
    single_dir=${coadd_dir/coadd/single}
    echo $single_dir
    header_dir=${coadd_dir/coadd/headers}
    file=${coadd_dir}/*.cut.fits
    base=$(basename ${file} _r.MP9602.V2.0.0A.swarp.cut.fits)
    num_pattern=${base:5}
    num=${num_pattern//[-._]/}
    ln -s ${coadd_dir}/${base}_r.MP9602.V2.0.0A.swarp.cut.fits ${output_dir}/tiles/CFIS_$num.fits
    ln -s ${coadd_dir}/${base}_r.MP9602.V2.0.0A.swarp.cut.weight.fits ${output_dir}/tiles/CFIS_$num.weight.fits
    ln -s ${coadd_dir}/${base}_r.MP9602.V2.0.0A.swarp.cut.flag.fits ${output_dir}/tiles/CFIS_$num.flag.fits
    # create exposure list, link files but exclude broken links
     for file in $(find ${single_dir}/*.sub.fits -type l -not -xtype l); do
        base=$(basename $file C.sub.fits)
        weight_link=${single_dir}/${base}C.weight.fits.gz
        if [ -L ${weight_link} ] ; then
            if [ -e ${weight_link} ] ; then
                echo $(basename $file C.sub.fits)
                ln -s ${single_dir}/${base}C.sub.fits ${output_dir}/exposures
                ln -s ${single_dir}/${base}C.weight.fits.gz ${output_dir}/exposures
                # remove header line 2 which has french accents that cause problems with astropy, change wcs keys
                awk 'NR!~/^(2)$/ {sub(/TAN/,"TPV"); print}' ${header_dir}/${base}.head > ${output_dir}/exposures/headers/${base}.head
                
            fi
         fi 
         done >> ${output_dir}/tiles/exposure_lists/exp_numbers-$num.txt
   
     #list=$(find ${single_dir}/*.sub.fits | awk -F/ '{print substr($NF,1,8)}'|uniq)
    #echo "${list}">exposure_lists/$num.txt
     
    done

   #uncompress step 
   gunzip -f ${output_dir}/exposures/*.gz
