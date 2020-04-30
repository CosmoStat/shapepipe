#!/bin/bash

# Description: Process one tile and all contributing exposures
#              on canfar
# Author: Martin Kilbinger <martin.kilbinger@cea.fr>
# Date: 03/2020
# Package: ShapePipe


### Set-up ###

## Variables

# Tile numbers

if [ $# == 0 ]; then
  echo "Usages:"
  echo "  bash sp_tile.bash TILE_ID_1 [TILE_ID_2 [...]]"
  echo "    TILE_ID = xxx.yyy"
  echo "    Example: sp_tile.bash 244.252"
  echo "  . sp_tile.bash -e"
  echo "    Assign environment variables"
  exit 1
fi

# Copy command line arguments
TILE_ARR=($@)

# For tar archives. Should be unique to each job
#ID=${TILE_ARR[0]}
ID=`echo ${TILE_ARR[@]} | tr ' ' '_'`

# For testing only use one exposure
ONE_EXP=0

## Paths

# VM home (do not modify)
export VM_HOME=/home/ubuntu

# Path to ShapePipe installation
export SP_ROOT=$VM_HOME/ShapePipe

## Path variables used in ShapePipe config files

# Run path and location of input image directories
export SP_RUN=`pwd`

# Config file path
export SP_CONFIG=$SP_RUN/GOLD

## Input and output paths used in config file

# Input
INPUT_TILES=$SP_RUN/input_tiles
INPUT_EXP=$SP_RUN/input_exposures

# Output
OUTPUT=$SP_RUN/output
OUTPUT_HEADERS=$SP_RUN/output_headers

# For tar archives
output_rel=`realpath --relative-to=. $OUTPUT`

# Temporary processing paths (e.g. for downloading images)
DOWNLOAD_EXP=$SP_RUN/download_exposures


## Other variables

# Stop on error
STOP=1

# Verbose mode (1: verbose, 0: quiet)
VERBOSE=1

if [ $VERBOSE == 1 ]; then
   vflag="-v"
else
   vflag=""
fi

# VCP options
# VCP without "vflag" to avoid output to stderr
export VCP_QUICK=1

if [ $VCP_QUICK == 1 ]; then
   qflag="--quick"
else
  qflag=""
fi
export CERTFILE=$VM_HOME/.ssl/cadcproxy.pem
export VCP="vcp $qflag --certfile=$CERTFILE"


## Functions

# Print string, executes command, and prints return value.
function command () {
   cmd=$1
   str=$2
   verbose=$3

   #RED='\033[0;31m'
   #GREEN='\033[0;32m'
   #NC='\033[0m' # No Color
   # Color escape characters show up in log files
   RED=''
   GREEN=''
   NC=''


   if [ $# == 3 ]; then
      if [ $verbose == 1 ]; then
           echo "$str: running '$cmd'"
      fi
      $cmd
   else
      if [ $verbose == 1 ]; then
         echo "$str: running '$cmd $4 \"$5 $6\"'"
      fi
      $cmd $4 "$5 $6"
   fi	
   res=$?

   if [ $verbose == 1 ]; then
      if [ $res == 0 ]; then
         echo -e "${GREEN}success, return value = $res${NC}"
      else
         echo -e "${RED}error, return value = $res${NC}"
         if [ $STOP == 1 ]; then
            echo "${RED}exiting 'sp_tile.bash', error in command '$cmd'${NC}"
            exit $res
	 else
            echo "${RED}continuing 'sp_tile.bash', error in command '$cmd'${NC}"
         fi
      fi
   fi

   return $res
}


# Run command. If error occurs, upload sp log files before stopping script.
command_sp() {
   cmd=$1
   str=$2
   verbose=$3
   id=$4

   STOP=0
   command "$1" "$2" "$3"
   res=$?
   if [ $res != 0 ]; then
      upload_logs $id $verbose
      echo "exiting 'sp_tile.bash', error in command '$cmd', log files for id=$id uploaded"
      exit $res
   fi

}

# Tar and upload files to vos
function upload() {
   base=$1
   shift
   ID=$1
   shift
   verbose=$1
   shift
   upl=("$@")

   n_upl=(`ls -l ${upl[@]} | wc`)
   if [ $n_upl == 0 ]; then
      if [ $STOP == 1 ]; then
         echo "Exiting script, no file found for '$base' tar ball"
         exit 3
      fi
   fi
   tar czf ${base}_${ID}.tgz ${upl[@]}
   command "$VCP ${base}_${ID}.tgz vos:cfis/cosmostat/kilbinger/results" "Upload $base results, $n_upl files in tar ball" "$verbose"
}

# Upload log files
function upload_logs() {
   id=$1
   verbose=$2

   upl="$output_rel/*/*/logs $output_rel/*/logs"
   upload "logs" "$id" "$verbose" "${upl[@]}"
}

# Print script variables
function print_env() {
   echo "*** Setting ***"
   echo "Data:"
   echo " TILE_ARR=${TILE_ARR[@]}"
   echo " ONE_EXP=$ONE_EXP"
   echo "Paths:"
   echo " VM_HOME=$VM_HOME"
   echo " SP_ROOT=$SP_ROOT"
   echo " SP_RUN=$SP_RUN"
   echo " INPUT_TILES=$INPUT_TILES"
   echo " INPUT_EXP=$INPUT_EXP"
   echo " OUTPUT=$OUTPUT"
   echo " OUTPUT_HEADERS=$OUTPUT_HEADERS"
   echo " DOWNLOAD_EXP=$DOWNLOAD_EXP"
   echo " SP_CONFIG=$SP_CONFIG"
   echo "Other variables:"
   echo " VCP=$VCP"
   echo " CERTFILE=$CERTFILE"
   echo " STOP=$STOP"
   echo " verbose=$VERBOSE"
   echo " LD_LIBRARY_PATH=$LD_LIBRARY_PATH"
   echo "***"
}


### Start ###

# Activate conda environment
echo "Activate conda 'shapepipe' environment"
source $VM_HOME/miniconda3/bin/activate shapepipe

print_env

# Extra stuff for canfar
export LD_LIBRARY_PATH=$VM_HOME/miniconda3/envs/shapepipe/lib

if [ "$1" == "-e" ]; then
   echo "Exiting"
   return 0
fi

echo "Start"

n_tile=${#TILE_ARR[@]}
echo "Processing $n_tile tile(s)"

## Create input and output directories

echo "Create directories for processing"
mkdir -p $SP_RUN
cd $SP_RUN
mkdir -p $DOWNLOAD_EXP
mkdir -p $INPUT_EXP
mkdir -p $INPUT_TILES
mkdir -p $OUTPUT
mkdir -p $OUTPUT_HEADERS


## Download and prepare input images


for TILE in ${TILE_ARR[@]}; do

   echo "Tile $TILE"

   # For cfis_field_select, remove '.' in file base name
   export SP_TILE=`echo $TILE | tr . -`

   # Get input tile and weight
   command "$VCP vos:cfis/tiles_DR2/CFIS.${TILE}.r.fits $INPUT_TILES/CFIS.${TILE}.r.fits" \
	   "Download tile image" "$VERBOSE"
   command "$VCP vos:cfis/tiles_DR2/CFIS.${TILE}.r.weight.fits.fz $INPUT_TILES/CFIS_weight-$SP_TILE.fits.fz" \
	   "Copy weight image from vos" "$VERBOSE"
   command "$SP_ROOT/scripts/python/cfis_weight_format.py -i $INPUT_TILES/CFIS_weight-$SP_TILE.fits.fz -o $INPUT_TILES/CFIS_weight-$SP_TILE.fits" \
	   "Transform weight" "$VERBOSE"

done


# Select exposures
command "$SP_ROOT/scripts/python/cfis_field_select.py -i $INPUT_TILES --tile -v -t exposure -o exp" \
	"Select exposures" "$VERBOSE"


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

### Debug: list tile files
ls -rtl $INPUT_TILES

# Download single exposure images, weights, flags
command "$SP_ROOT/scripts/python/cfis_download_images.py -i exp.txt -o $DOWNLOAD_EXP $vflag -t exposure --in_number_only -q --certfile $CERTFILE" \
	"Download exposure images" "$VERBOSE"
command "$SP_ROOT/scripts/python/cfis_download_images.py -i exp.txt -o $DOWNLOAD_EXP $vflag -t exposure_weight.fz --in_number_only --certfile $CERTFILE" \
	"Download exposure weights" "$VERBOSE"
command "$SP_ROOT/scripts/python/cfis_download_images.py -i exp.txt -o $DOWNLOAD_EXP $vflag -t exposure_flag.fz --in_number_only --certfile $CERTFILE" \
	"Download flags" "$VERBOSE"

### Debug: Check why sometimes files are missing even though cfis_download_images.py exists without error
ls -rtl $DOWNLOAD_EXP/*.fits.fz
n_downl=(`ls -l $DOWNLOAD_EXP/*.fits.fz | wc`)
n_downl_X=`perl -e 'print '$n_exp' * 3'`
echo "$n_downl images downloaded, expected $n_downl_X"
if [ $n_downl != $n_downl_X ]; then
   if [ $STOP == 1 ]; then
      echo "Exiting script, number of downloaded files does not match"
      exit 98
   fi
fi

# Set links to exposures, weights, and flags
for i in `cat $SP_RUN/exp.txt`; do
   src=$DOWNLOAD_EXP/${i}p.fits.fz
   trg=$INPUT_EXP/image-$i.fitsfz
   cmd="ln -sf $src $trg"
   $cmd

   src=$DOWNLOAD_EXP/${i}p.weight.fits.fz
   trg=$INPUT_EXP/weight-$i.fitsfz
   cmd="ln -sf $src $trg"
   $cmd

   src=$DOWNLOAD_EXP/${i}p.flag.fits.fz
   trg=$INPUT_EXP/flag-$i.fitsfz
   cmd="ln -sf $src $trg"
   $cmd
done

# Download config files
$VCP vos:cfis/cosmostat/kilbinger/GOLD .


### Run pipeline

## Exposures

# Run all modules: Leads to unresolved thread error
#$SP_ROOT/shapepipe_run.py -c $SP_CONFIG/config_exp.ini

# Split up, merge headers, and mask
command_sp "$SP_ROOT/shapepipe_run.py -c $SP_CONFIG/config_exp_SpMeMa.ini" "Run ShapePipe 1/5 (exp: split, merge headers, mask)" "$VERBOSE" "$ID"

# Source-extract and select stars
command_sp "$SP_ROOT/shapepipe_run.py -c $SP_CONFIG/config_exp_SxSt.ini" "Run ShapePipe 2/5 (exp: extract, select)" "$VERBOSE" "$ID"

# Create PSF model
command_sp "$SP_ROOT/shapepipe_run.py -c $SP_CONFIG/config_exp_Ps.ini" "Run ShapePipe 3/5 (exp: create PSF)" "$VERBOSE" "$ID"

# Interpolate PSF model for validation
command_sp "$SP_ROOT/shapepipe_run.py -c $SP_CONFIG/config_exp_Pi.ini" "Run ShapePipe 4/5 (exp: interpolate PSF)" "$VERBOSE" "$ID"


# The following are very a bad hacks to get additional input files
input_psfex=`find . -name star_selection-*.psf | head -n 1`
command_sp "ln -s `dirname $input_psfex` input_psfex" "Link psfex output" "$VERBOSE" "$ID"

input_split_exp=`find output -name flag-*.fits | head -n 1`
command_sp "ln -s `dirname $input_split_exp` input_split_exp" "Link split_exp output" "$VERBOSE" "$ID"

input_sextractor=`find . -name sexcat_sexcat-*.fits | head -n 1`
command_sp "ln -s `dirname $input_sextractor` input_sextractor" "Link sextractor output" "$VERBOSE" "$ID"

## Tiles
command_sp "$SP_ROOT/shapepipe_run.py -c $SP_CONFIG/config_tile.ini" "Run ShapePipe 5/5 (tile: everything)" "$VERBOSE" "$ID"


# module and pipeline log files
upload_logs "$ID" "$VERBOSE"


## Upload results

# psfex for diagnostics, validation with leakage
upl=$output_rel/*/psfex_runner/output/star_selection-*
upload "psfex" "$ID" "$VERBOSE" "${upl[@]}"

# psefxinterp for validation with residuals, rho stats
upl=$output_rel/*/psfexinterp_runner/output/validation_psf*
upload "psfexinterp" "$ID" "$VERBOSE" "${upl[@]}"

# SETools stats and plots
upl=$output_rel/*/setools_runner/output/stat/*
upload "setools_stat" "$ID" "$VERBOSE" "${upl[@]}"
upl=$output_rel/*/setools_runner/output/plot/*
upload "setools_plot" "$ID" "$VERBOSE" "${upl[@]}"

# shapes
if [ -e $output_rel/*/ngmix_runner/output ]; then
   n_file=(`ls -l $output_rel/*/ngmix_runner/output | wc`)
   if [ "$n_file" == 1 ]; then
      if [ $STOP == 1 ]; then
         echo "Existing script, no ngmix FITS file found in ngmix output dir"
         exit 97
      fi
   else
      upl=$output_rel/*/ngmix_runner/output/ngmix-*
      n_upl=(`ls -l ${upl[@]} | wc`)
      if [ $n_upl != $n_tile ]; then
         echo  "Warning: number of ngmix files found ($n_upl) does not match expectation ($n_tile)"
      fi
      upload "ngmix" "$ID" "$VERBOSE" "${upl[@]}"
   fi
else
   if [ $STOP == 1 ]; then
      echo "Existing script, no ngmix output dir found"
      exit 98
   fi
fi

echo "End"
