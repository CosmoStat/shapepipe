#!/usr/bin/bash

# Usage
# ~/astro/repositories/github/shapepipe/scripts/sh/curl_canfar.sh 0.9 shapepipe/scripts/sh/init_run_exlusive_canfar.sh ID ind

SSL=~/.ssl/cadcproxy.pem
N_SMP=2
SESSION=https://ws-uv.canfar.net/skaha/v0/session
RESOURCES="ram=4&cores=$N_SMP"
IMAGE=images.canfar.net/unions/shapepipe
NAME=shapepipe

# version of image on canfar, e.g. 0:7, 0:8
version=$1

# command on canfar, e.g. shapepipe/scripts/sh/init_run_exclusive_canfar.sh
cmd=$2

# Kind ("tile" or "exp")
kind=$3

# Image ID; has to be last argument to work with xargs
ID=$4

# command line arguments for remote script:
# collect into string
arg="-j $JOB -e $ID -n $N_SMP -k $kind"

ID=`curl -E $SSL $SESSION?$RESOURCES -d "image=$IMAGE:$version" -d "name=${NAME}" -d "cmd=$cmd" --data-urlencode "args=$arg"`
echo $ID >> IDs.txt
