#!/usr/bin/env bash

# Name: canfar_sp.bash
# Description: Process one or more tiles with all
#              contributing exposures on canfar.
#              This is the job submission script for
#              the canfar batch system. Can also be
#              called in interactive mode on a virtual
#              machine.
# Author: Martin Kilbinger <martin.kilbinger@cea.fr>
# Date: 03/2020
# Package: shapepipe


### Set-up ###

## Variables

# Tile numbers

if [ $# == 0 ]; then
  echo "Usages:"
  echo "  bash canfar_sp.bash TILE_ID_1 [TILE_ID_2 [...]]"
  echo "    TILE_ID = xxx.yyy"
  echo "    Examples:"
  echo "      canfar_sp.bash 244.252"
  echo "      canfar_sp.bash 244.252 239.293"
  echo "  . canfar_sp.bash -e"
  echo "    Assign environment variables"
  exit 1
fi

# Copy command line arguments
TILE_ARR=($@)

# For tar archives. Should be unique to each job
export ID=`echo ${TILE_ARR[@]} | tr ' ' '_'`


## Paths

# VM home (do not modify)
export VM_HOME=/home/ubuntu

## Path variables used in shapepipe config files

# Run path and location of input image directories
export SP_RUN=`pwd`

# Config file path
export SP_CONFIG=$SP_RUN/cfis

# Input tile numbers ASCII file
export TILE_NUMBERS_PATH=tile_numbers.txt

# Output
OUTPUT=$SP_RUN/output

# For tar archives
output_rel=`realpath --relative-to=. $OUTPUT`


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
            echo "${RED}exiting 'canfar_sp.bash', error in command '$cmd'${NC}"
            exit $res
         else
            echo "${RED}continuing 'canfar_sp.bash', error in command '$cmd'${NC}"
         fi
      fi
   fi

   return $res
}

# Run shapepipe command. If error occurs, upload sp log files before stopping script.
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
      echo "exiting 'canfar_sp.bash', '$cmd' returned $res, log files for id=$id uploaded"
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
   echo "Paths:"
   echo " VM_HOME=$VM_HOME"
   echo " SP_RUN=$SP_RUN"
   echo " TILE_NUMBERS_PATH=$TILE_NUMBERS_PATH"
   echo " OUTPUT=$OUTPUT"
   echo " SP_CONFIG=$SP_CONFIG"
   echo "Other variables:"
   echo " VCP=$VCP"
   echo " CERTFILE=$CERTFILE"
   echo " qflag=$qflag"
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

# Create input and output directories
echo "Create directories for processing"
mkdir -p $SP_RUN
cd $SP_RUN
mkdir -p $OUTPUT

# Write tile numbers to ASCII input file
rm -rf $TILE_NUMBERS_PATH
for TILE in ${TILE_ARR[@]}; do
   echo $TILE >> $TILE_NUMBERS_PATH
done

# Download config files
$VCP vos:cfis/cosmostat/kilbinger/cfis .

### Run pipeline

## Prepare images

# Get tiles
command_sp "shapepipe_run -c $SP_CONFIG/config_get_tiles.ini" "Run shapepipe (prepare: get tiles)" "$VERBOSE" "$ID"

# Uncompress tile weights
command_sp "shapepipe_run -c $SP_CONFIG/config_unfz_w.ini" "Run shapepipe (prepare: uncompress tile weights)" "$VERBOSE" "$ID"

# Find exposures
command_sp "shapepipe_run -c $SP_CONFIG/config_find_exp.ini" "Run shapepipe (prepare: find exposures)" "$VERBOSE" "$ID"

# Get exposures
command_sp "shapepipe_run -c $SP_CONFIG/config_get_exp.ini" "Run shapepipe (prepare: get exposures)" "$VERBOSE" "$ID"

## Exposures

# Run all modules
command_sp "shapepipe_run -c $SP_CONFIG/config_exp.ini" "Run shapepipe (exp)" "$VERBOSE" "$ID"


# The following are very a bad hacks to get additional input files
input_psfex=`find . -name star_split_ratio_80-*.psf | head -n 1`
command_sp "ln -s `dirname $input_psfex` input_psfex" "Link psfex output" "$VERBOSE" "$ID"

input_split_exp=`find output -name flag-*.fits | head -n 1`
command_sp "ln -s `dirname $input_split_exp` input_split_exp" "Link split_exp output" "$VERBOSE" "$ID"

input_sextractor=`find . -name sexcat_sexcat-*.fits | head -n 1`
command_sp "ln -s `dirname $input_sextractor` input_sextractor" "Link sextractor output" "$VERBOSE" "$ID"


## Tiles

# Everything up to shapes
command_sp "shapepipe_run -c $SP_CONFIG/config_tile_MaSxPsViSmVi.ini" "Run shapepipe (tile: up to ngmix)" "$VERBOSE" "$ID"

# Shapes, run 8 parallel processes
for k in $(seq 1 8); do
    command_sp "shapepipe_run -c $SP_CONFIG/config_tile_Ng${k}u.ini" "Run shapepipe (tile: ngmix $k)" "$VERBOSE" "$ID" &
done
wait

# Merge separated shapes catalogues
command_sp "shapepipe_run -c $SP_CONFIG/config_merge_sep_cats.ini" "Run shapepipe (tile: merge sep cats)" "$VERBOSE" "$ID"

# Create final shape catalogue by merging all tile information
command_sp "shapepipe_run -c $SP_CONFIG/config_make_cat.ini" "Run shapepipe (tile: create final cat)" "$VERBOSE" "$ID"


## Upload results

# module and pipeline log files
upload_logs "$ID" "$VERBOSE"

# psfex for diagnostics, validation with leakage
# psefxinterp for validation with residuals, rho stats
# SETools masks (selection), stats and plots
# Final shape catalog

NAMES=(
        "psfex"
        "psfexinterp_exp"
        "setools_mask"
        "setools_stat"
        "setools_plot"
        "final_cat"
     )
DIRS=(
        "*/psfex_runner/output"
        "*/psfexinterp_runner/output"
        "*/setools_runner/output/mask"
        "*/setools_runner/output/stat"
        "*/setools_runner/output/plot"
        "*/make_catalog_runner/output"
     )
PATTERNS=(
        "star_split_ratio_80-*"
        "validation_psf*"
        "*"
        "*"
        "*"
        "final_cat-*"
        )

for n in "${!NAMES[@]}"; do
    name=${NAMES[$n]}
    dir=${DIRS[$n]}
    pattern=${PATTERNS[$n]}
    upl=$output_rel/$dir/$pattern
    upload "$name" "$ID" "$VERBOSE" "${upl[@]}"
done

echo "End"
