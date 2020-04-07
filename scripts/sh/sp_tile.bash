#!/bin/bash

# Description: Process one tile and all contributing exposures
# Author: Martin Kilbinger <martin.kilbinger@cea.fr>
# Date: 03/2020
# Package: ShapePipe


### Set-up ###

## Variables

# Tile numbers

if [ $# == 0 ]; then
  echo "Usage: sp_tile.bash TILE_ID_1 [TILE_ID_2 [...]]"
  echo "  TILE_ID = xxx.yyy"
  echo "  Example: sp_tile.bash 244.252"
  exit 1

  # Or default (test) tile
  #TILE_ARR=(619.241)
fi

# Copy command line arguments
TILE_ARR=($@)


# For testing only use one exposure
ONE_EXP=0

## VM-specific paths, do not modify
# VM home
export VM_HOME=/home/ubuntu

# VM path to ShapePipe installation, used in config files
export SP_ROOT=$VM_HOME/ShapePipe
## End VM-specific paths

# Processing paths (used in the pipeline config files)

# Run path, used in config files.
export SP_RUN=`pwd`

# Input
export INPUT_TILES=$SP_RUN/input_tiles
export INPUT_EXP=$SP_RUN/input_exposures

# Output
export OUTPUT=$SP_RUN/output
export OUTPUT_HEADERS=$SP_RUN/output_headers
#export PSF_VALIDATION=$SP_RUN/psf_validation

# Temporary processing paths (e.g. for downloading images)
DOWNLOAD_EXP=$SP_RUN/download_exposures

# Config file path
export SP_CONFIG=$SP_RUN/GOLD

## Other variables

# Verbose mode (1: verbose, 0: quiet)
verbose=1

if [ $verbose == 1 ]; then
   vflag="-v"
else
   vflag=""
fi

# Command shortcuts
# VCP without "vflag" to avoid to output to stderr
VCP="vcp --quick --certfile=$VM_HOME/.ssl/cadcproxy.pem"


## Functions

# Print string, executes command, and prints return value.
# Printing if verbose=1
function command () {
   cmd=$1
   str=$2

   #RED='\033[0;31m'
   #GREEN='\033[0;32m'
   #NC='\033[0m' # No Color
   # Escape characters show up in log files
   RED=''
   GREEN=''
   NC=''


   if [ $# == 2 ]; then
   	[ $verbose ] &&  echo "$str: running '$cmd'"
   	$cmd
   else
   	[ $verbose ] &&  echo "$str: running '$cmd $3 \"$4 $5\"'"
	$cmd $3 "$4 $5"
   fi	
   res=$?

   if [ $verbose ]; then
      if [ $res == 0 ]; then
         echo -e "${GREEN}success, return value = $res${NC}"
      else
         echo -e "${RED}error, return value = $res${NC}"
      fi
   fi
}

# Print script variables
function print_env() {
   echo "*** Setting ***"
   echo "Data:"
   echo " TILE_ARR = ${TILE_ARR[@]}"
   echo " ONE_EXP = $ONE_EXP"
   echo "Paths:"
   echo " VM_HOME = $VM_HOME"
   echo " SP_ROOT = $SP_ROOT"
   echo " SP_RUN = $SP_RUN"
   echo " INPUT_TILES = $INPUT_TILES"
   echo " INPUT_EXP = $INPUT_EXP"
   echo " OUTPUT = $OUTPUT"
   echo " OUTPUT_HEADERS = $OUTPUT_HEADERS"
   #echo " PSF_VALIDATION = $PSF_VALIDATION"
   echo " DOWNLOAD_EXP = $DOWNLOAD_EXP"
   echo " SP_CONFIG = $SP_CONFIG"
   echo "Other variables:"
   echo " VCP = $VCP"
   echo " verbose = $verbose"
   echo " LD_LIBRARY_PATH = $LD_LIBRARY_PATH"
   echo "***"
}


### Start ###

echo "Start"


echo "Processing ${#TILE_ARR[@]} tile(s)"

# Activate conda environment
echo "Activate conda 'shapepipe' environment"
source $VM_HOME/miniconda3/bin/activate shapepipe

# Extra stuff for canfar
export LD_LIBRARY_PATH=$VM_HOME/miniconda3/envs/shapepipe/lib

print_env

## Create input and output directories

echo "Create directories for processing"
mkdir -p $SP_RUN
cd $SP_RUN
mkdir -p $DOWNLOAD_EXP
mkdir -p $INPUT_EXP
mkdir -p $INPUT_TILES
mkdir -p $OUTPUT
mkdir -p $OUTPUT_HEADERS
#mkdir -p $PSF_VALIDATION


## Download and prepare input images


for TILE in ${TILE_ARR[@]}; do

   echo "Tile $TILE"

   # For cfis_field_select, remove '.' in file base name
   export SP_TILE=`echo $TILE | tr . -`

   # Get input tile and weight
   command "$VCP vos:cfis/tiles_DR2/CFIS.${TILE}.r.fits $INPUT_TILES/CFIS.${TILE}.r.fits" \
	   "Download tile image"
   command "$VCP vos:cfis/tiles_DR2/CFIS.${TILE}.r.weight.fits.fz $INPUT_TILES/CFIS_weight-$SP_TILE.fits.fz" \
	   "Copy weight image from vos"
   command "$SP_ROOT/scripts/python/cfis_weight_format.py -i $INPUT_TILES/CFIS_weight-$SP_TILE.fits.fz -o $INPUT_TILES/CFIS_weight-$SP_TILE.fits" \
	   "Transform weight"

done


# Select exposures
command "$SP_ROOT/scripts/python/cfis_field_select.py -i input_tiles --tile -v -t exposure -o exp" \
	"Select exposures"


### Only leave one line in exposure list
if [ $ONE_EXP == 1 ]; then
	echo "=== For testing: Use only one exposure! ==="
	head -n 1 exp.txt > exp1.txt
	mv exp1.txt exp.txt
fi


# Count number of exposures to process
n_exp=(`wc exp.txt`)
echo "$n_exp exposures to be processed"

# Rename tile, needs to be done after previous 'cfis_field_select.py' call
mv $INPUT_TILES/CFIS.${TILE}.r.fits $INPUT_TILES/CFIS_image-$SP_TILE.fits

# Download single exposure images, weights, flags
command "$SP_ROOT/scripts/python/cfis_download_images.py -i exp.txt -o $DOWNLOAD_EXP $vflag -t exposure --in_number_only" \
	"Download exposure images" "--vcp" "$VCP"
command "$SP_ROOT/scripts/python/cfis_download_images.py -i exp.txt -o $DOWNLOAD_EXP $vflag -t exposure_weight.fz --in_number_only" \
	"Download exposure weights" "--vcp" "$VCP"
command "$SP_ROOT/scripts/python/cfis_download_images.py -i exp.txt -o $DOWNLOAD_EXP $vflag -t exposure_flag.fz --in_number_only" \
	"Download flags" "--vcp" "$VCP"

n_downl=(`ls -l $DOWNLOAD_EXP/*.fits.fz | wc`)
n_downl_X=`perl -e 'print '$n_exp' * 3'`
echo "$n_downl images downloaded, expected $n_downl_X"

# Set links to exposures, weights, and flags
for i in `cat $SP_RUN/exp.txt`; do
   src=$DOWNLOAD_EXP/${i}p.fits.fz
   echo $src
   trg=$INPUT_EXP/image-$i.fitsfz
   cmd="ln -sf $src $trg"
   #command "$cmd" "symlink $src <- $trg"
   $cmd

   src=$DOWNLOAD_EXP/${i}p.weight.fits.fz
   trg=$INPUT_EXP/weight-$i.fitsfz
   cmd="ln -sf $src $trg"
   #command "$cmd" "symlink $src <- $trg"
   $cmd

   src=$DOWNLOAD_EXP/${i}p.flag.fits.fz
   trg=$INPUT_EXP/flag-$i.fitsfz
   cmd="ln -sf $src $trg"
   #command "$cmd" "$symlink $src <- $trg"
   $cmd
done

# Download config files
$VCP vos:cfis/cosmostat/kilbinger/GOLD .


### Run pipeline

## Exposures

$SP_ROOT/shapepipe_run.py -c $SP_CONFIG/config_exp.ini


## Upload results

# psfex
upl=$OUTPUT/*/psfex_runner/output/star_selection-*
n_upl=(`ls -l $upl | wc`)
command "$VCP $upl vos:cfis/cosmostat/kilbinger/psfex" "Upload psfex results, $n_up files"

# psefxinterp for validation
upl=$OUTPUT/*/psfexinterp_runner/output/validation_psf*
n_upl=(`ls -l $upl | wc`)
command "$VCP $upl vos:cfis/cosmostat/kilbinger/psfexinterp" "Upload psfexinterp results, $n_up files"

# module log files
upl=$OUTPUT/*/*/logs
n_upl=(`ls -l $upl | wc`)
command "$VCP $upl vos:cfis/cosmostat/kilbinger/logs" "Upload module logs, $n_up files"

upl=$OUTPUT/*/logs
n_upl=(`ls -l $upl | wc`)
command "$VCP $upl vos:cfis/cosmostat/kilbinger/logs" "Upload shapepipe logs, $n_up files"
echo "End"
