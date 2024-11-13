#!/usr/bin/env bash

# Global variables
SSL=~/.ssl/cadcproxy.pem
SESSION=https://ws-uv.canfar.net/skaha/v0/session
IMAGE=images.canfar.net/unions/shapepipe
NAME=shapepipe

# :TODO: Not working
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $HOME/shapepipe/scripts/sh/functions.sh

# Command line arguments

## Default values
job=-1
psf="psfex"
ID=-1
file_IDs=-1
N_SMP=1
version="1.1"
cmd_remote="$HOME/shapepipe/scripts/sh/init_run_exclusive_canfar.sh"
batch_max=200
dry_run=0
mh_local=0
debug_out="-1"

pat="- "

## Help string
usage="Usage: $(basename "$0") -j JOB -[e ID |-f file_IDs] -k KIND [OPTIONS]
\n\nOptions:\n
   -h\tthis message\n
   -j, --job JOB\tRunning JOB, bit-coded\n
   -e, --exclusive ID
    \timage ID\n
   -f, --file_IDs path
    \tfile containing IDs\n
   -p, --psf MODEL\n
    \tPSF model, one in ['psfex'|'mccd'], default='$psf'\n
   -m, --mh_local MH\n
    \tmerged header file local (MH=0) or global (MH=1); default is $mh_local\n
   -N, --N_SMP N_SMOp\n
    \tnumber of jobs (SMP mode only), default=$N_SMP\n
   -V, --version\n
    \tversion of docker image, default='$version'\n
   -C, --command_remote\n
    \tremote command to run on canfar, default='$cmd_remote'\n
   -b, --batch_max\n
    \tmaximum batch size = number of jobs run simultaneously, default=$batch_max\n
   --debug_out PATH\n
    \tdebug output file PATH, default not used\n
   -n, --dry_run LEVEL\n
    \tdry run, from LEVEL=2 (no processing) to 0 (full run)\n
"

## Help if no arguments
if [ -z $1 ]; then
  echo -ne $usage
  exit 1
fi

## Parse command line
while [ $# -gt 0 ]; do
  case "$1" in
    -h)
      echo -ne $usage
      exit 0
      ;;
    -j|--job)
      job="$2"
      shift
      ;;
    -p|--psf)
      psf="$2"
      shift
      ;;
    -m|--mh_local)
      mh_local="$2"
      shift
      ;;
    -e|--exclusive)
      ID="$2"
      shift
      ;;
    -f|--file_IDs)
      file_IDs="$2"
      shift
      ;;
    -N|--N_SMP)
      N_SMP="$2"
      shift
      ;;
    -b|--batch_max)
      batch_max="$2"
      shift
      ;;
    --debug_out)
      debug_out="$2"
      shift
      ;;
    -n|--dry_run)
      dry_run="$2"
      shift
      ;;
  esac
  shift
done

## Check options                                                                 
if [ "$job" == "-1" ]; then                                                     
  echo "No job indicated, use option -j"                                        
  exit 2                                                                        
fi                                                                              
                                                                                
if [ "$ID" == "-1" ] && [ "$file_IDs" == "-1" ]; then                                               
  echo "No image ID(s) indicated, use option -e ID or -f file_IDs"                                   
  exit 3                                                                        
fi                                                                              

if [ "$psf" != "psfex" ] && [ "$psf" != "mccd" ]; then
  echo "PSF (option -p) needs to be 'psfex' or 'mccd'"
  exit 4
fi
                                                                                
if [ "$dry_run" != 0 ] && [ "$dry_run" != 1 ] && [ "$dry_run" != 2 ]; then
  echo "Invalid dry_run option, allowed are 0, 1, and 2"
  exit 5
fi

if [ "$debug_out" != "-1" ]; then
  echo "${pat}Starting $(basename "$0")" >> $debug_out
  echo $pat`date`$ >> $debug_out
fi

. /opt/conda/etc/profile.d/conda.sh
conda activate shapepipe
if [ "$debug_out"  != "-1" ]; then
    echo "${pat}conda prefix = ${CONDA_PREFIX}" >> $debug_out
fi

# command line arguments for remote script:
# collect into string

if [ "$dry_run" == "1" ]; then
  arg_dry_run="-n $dry_run"
else
  arg_dry_run=""
fi

RESOURCES="ram=4&cores=$N_SMP"
dir=`pwd`


function submit_batch() {
  path=$1

  # Paralle call of curl to speed up submission
  #n_para=8
  #cat $path | xargs -I {} -P $n_para bash -c '
    #source $HOME/shapepipe/scripts/sh/functions.sh
    #ID="{}"
    #IDt=$(echo $ID | tr "." "-")
    #my_name="SP-${patch}-J${job}-${IDt}"
    #call_curl $my_name $dry_run $debug_out
  #'

  for ID in `cat $path`; do
    IDt=`echo $ID | tr "." "-"`
    my_name="SP-${patch}-J${job}-${IDt}"
    call_curl $my_name $dry_run $debug_out
  done
}

batch=50
sleep=75

((n_thresh=batch_max-batch))


if [ "$dry_run" == 2 ]; then

  # Do not call curl (dry run = 2)
  echo "Running command dry run:"

  if [ "$ID" == "-1" ]; then


    # Submit file (dry run = 2)
    for ID in `cat $file_IDs`; do
      IDt=`echo $ID | tr "." "-"`
      my_name="SP-${patch}-J${job}-${IDt}"
      call_curl $my_name $dry_run $debug_out
    done

  else

    # Submit image (dry run = 2)
    IDt=`echo $ID | tr "." "-"`
    my_name="SP-${patch}-J${job}-${IDt}"
    call_curl $my_name $dry_run $debug_out

  fi

else

  # Call curl
  rm -rf session_IDs.txt session_image_IDs.txt

  if [ "$ID" == "-1" ]; then

    # Submit file
    n_jobs=`cat $file_IDs | wc -l`
    if [ "$n_jobs" -gt "$batch_max" ]; then

      # Split into batches 
      prefix="${file_IDs}_split_"
      split -d -l $batch $file_IDs $prefix
      n_split=`ls -l $prefix* | wc -l`
      echo "Split '$file_IDs' into $n_split batches of size $batch"

      count=1
      n_running=`stats_jobs_canfar.sh`
      for batch in $prefix*; do
        echo "Number of running jobs = $n_running"
        echo "Submitting batch $batch ($count/$n_split)"
        echo -ne "\033]0;curl patch=$patch job=$job $count/$n_split\007"
        submit_batch $batch
        ((count=count+1))

        n_running=`stats_jobs_canfar.sh`

        while [ "$n_running" -gt "$n_thresh" ]; do
          echo "Wait for #jobs = $n_running jobs to go < $n_thresh ..."
          sleep $sleep
          n_running=`stats_jobs_canfar.sh`
        done

      done

    else

      # Submit entire file (single batch)
      echo "Submit '$file_IDs' in single batch"
      submit_batch $file_IDs

    fi

  else

    # Submit image
    IDt=`echo $ID | tr "." "-"`
    my_name="SP-${patch}-J${job}-${IDt}"
    call_curl $my_name $dry_run $debug_out

  fi

fi

echo "Done $(basename "$0")" 

if [ "$debug_out" != "-1" ]; then
  echo "${pat}End $(basename "$0")" >> $debug_out
fi
