#!/usr/bin/env bash

# Name: submit_selection.sh
# Author: Martin Kilbinger (2020), martin.kilbinger@cea.fr
# Description: Submits jobs to canfar, each jobs processes
#  single exposures for one tile
# Package: ShapePipe

## Variables

# Job file for one tile
SP_ROOT=$HOME/shapepipe
sp_job="$SP_ROOT/scripts/sh/canfar_sp.bash"

# Image/VM ID
image="ShapePipe2-mk-20200619"


if [ "$1" == "-h" ]; then
  echo "Usage: submit_selection.sh [-n]"
  echo "Options:"
  echo "  -n       dry run"
  exit 0
fi


if [ "$1" == "-n" ]; then
  dry_run=1
  dry_str=" (dry run)"
else
  dry_run=0
  dry_str=""
fi


# Functions

function add_queue() {
   local job_file=$1

   echo "queue" >> $job_file
   echo >> $job_file
}

function create_job_file() {
   local tile_ID=$1
   local job_file=$2

   echo "executable     = $sp_job" > $job_file
   echo >> $job_file

   add_to_job_file $tile_ID $job_file

   finalize_job_file $job_file
   add_queue $job_file
}


function add_to_job_file () {
   local tile_ID=$1
   local job_file=$2

   echo "arguments      = $tile_ID" >> $job_file
   echo "output         = log_canfar_sp_$tile_ID.out" >> $job_file
   echo "error          = log_canfar_sp_$tile_ID.err" >> $job_file
   echo "log            = log_canfar_sp_$tile_ID.log" >> $job_file
}


function finalize_job_file() {
   local job_file=$1

   echo "request_cpus   = 8" >> $job_file

   # Cannot be larger than VM available RAM
   echo "request_memory = 19G" >> $job_file
   echo "request_disk   = 100G" >> $job_file

   echo >> $job_file
}

if [ $dry_run == 0 ]; then
   # TODO: Check whether the following directories exist
   # (e.g. error msg of vmkdir), remove previous run
   # results
   vmkdir vos:cfis/cosmostat/kilbinger/results
fi


# Input tile list file
tile_list_file="../w3_tile.txt"

# Output tile ID list file
tile_ID_list="tile_id_list.txt"

# Create tile ID list
cat $tile_list_file | perl -ane '/CFIS\.(\d{3}\.\d{3})/; print "$1\n"' > $tile_ID_list


# Create job file
job_file="job_tile.sh"
echo "executable     = $sp_job" > $job_file
echo >> $job_file

echo "output         = log_canfar_sp_\$(arguments).out" >> $job_file
echo "error          = log_canfar_sp_\$(arguments).err" >> $job_file
echo "log            = log_canfar_sp_\$(arguments).log" >> $job_file
echo >> $job_file
finalize_job_file $job_file
echo >> $job_file
echo "queue arguments from (" >> $job_file

# Go through tile IDs and add job
while read -r tile_ID; do
   echo "$tile_ID" >> $job_file
done < $tile_ID_list
 echo ")" >> $job_file

# Submit
cmd="canfar_submit $job_file $image c8-90gb-186"
echo "Running $cmd$dry_str"
if [ $dry_run == 0 ]; then
   $cmd
   condor_q
fi
