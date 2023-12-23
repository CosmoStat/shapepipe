#!/usr/bin/env bash

# Usage
# ~/astro/repositories/github/shapepipe/scripts/sh/curl_canfar.sh 0.9 shapepipe/scripts/sh/init_run_exlusive_canfar.sh ID kind job

# Command line arguments
## Default values
SSL=~/.ssl/cadcproxy.pem
SESSION=https://ws-uv.canfar.net/skaha/v0/session
IMAGE=images.canfar.net/unions/shapepipe
NAME=shapepipe

job=-1
ID=-1
N_SMP=1
kind=-1
version=0.9
cmd_remote="shapepipe/scripts/sh/init_run_exclusive_canfar.sh"
dry_run=0

# TODO psf

## Help string
usage="Usage: $(basename "$0") -j JOB -e ID  -k KIND [OPTIONS]
\n\nOptions:\n
   -h\tthis message\n
   -j, --job JOB\tRUnning JOB, bit-coded\n
   -e, --exclusive ID
    \timage ID\n
   -p, --psf MODEL\n
    \tPSF model, one in ['psfex'|'mccd'], default='$psf'\n
   -k, --kind KIND\n
    \timage kind, allowed are 'tile' and 'exp'\n
   -N, --N_SMP N_SMOp\n
    \tnumber of jobs (SMP mode only), default from original config files\n
   -V, --version\n
    \tversion of docker image, default='$version'\n
   -C, --command_remote\n
    \tremote command to run on canfar, default='$cmd_remote'\n
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
    -e|--exclusive)
      ID="$2"
      shift
      ;;
    -N|--N_SMP)
      n_SMP="$2"
      shift
      ;;
    -k|--kind)
      kind="$2"
      shift
      ;;
    -n|--dry_run)
      dry_run="$2"
      shift
      ;;
  esac
  shift
done

# Check options                                                                 
if [ "$job" == "-1" ]; then                                                     
  echo "No job indicated, use option -j"                                        
  exit 2                                                                        
fi                                                                              
                                                                                
if [ "$exclusive" == "-1" ]; then                                               
  echo "No image ID indicated, use option -e"                                   
  exit 3                                                                        
fi                                                                              
                                                                                
if [ "kind" == "-1" ]; then                                                     
  echo "No image kind indicated, use option -k"                                 
  exit 4                                                                        
fi

if [ "$dry_run" != 0 ] && [ "$dry_run" != 1 ] && [ "$dry_run" != 2 ]; then
  echo "Invalid dry_run option, allowed are 0, 1, and 2"
  exit 5
fi

# command line arguments for remote script:
# collect into string

if [ "$dry_run" == "1" ]; then
  arg_dry_run="-n $dry_run"
else
  arg_dry_run=""
fi

RESOURCES="ram=4&cores=$N_SMP"
arg="-j $job -e $ID -N $N_SMP -k $kind $arg_dry_run"


if [ "$dry_run" == 2 ]; then
  echo "Running command dry run:"
  echo curl -E $SSL $SESSION?$RESOURCES -d \"image=$IMAGE:$version\" -d \"name=${NAME}\" -d \"cmd=$cmd_remote\" --data-urlencode \"args=$arg\"
else
  session=`curl -E $SSL $SESSION?$RESOURCES -d "image=$IMAGE:$version" -d "name=${NAME}" -d "cmd=$cmd_remote" --data-urlencode "args=$arg"`
  echo $session >> session_IDs.txt
  echo "$session $ID" >> session_image_IDs.txt
fi
