#!/usr/bin/bash

# Usage
# ~/astro/repositories/github/shapepipe/scripts/sh/curl_canfar.sh 0.8 shapepipe/scripts/sh/init_run_exlusive_canfar.sh ID exp NCORE

SSL=~/.ssl/cadcproxy.pem
NCORE=2
SESSION=https://ws-uv.canfar.net/skaha/v0/session
RESOURCES="ram=4&cores=$NCORE"
IMAGE=images.canfar.net/unions/shapepipe
NAME=shapepipe

# version of image on canfar, e.g. 0:7, 0:8
version=$1

# command on canfar, e.g. shapepipe/scripts/sh/init_run_exclusive_canfar.sh
cmd=$2

# Image ID
ID=$3

# command line arguments
arg="$ID $NCORE exp"

#arg_an=`echo "$arg" | tr '_' 'X' | tr '.' 'X'`

echo
echo "Start headless container"
echo "========================"
ID=`curl -E $SSL $SESSION?$RESOURCES -d "image=$IMAGE:$version" -d "name=${NAME}" -d "cmd=$cmd" --data-urlencode "args=$arg"`
echo $ID

echo
echo "Events (incl. errors)"
echo "====================="
cmd="curl -E $SSL $SESSION/$ID?view=events"
echo $cmd
$cmd

echo
echo "Logs (incl. stdout)"
echo "==================="
cmd="curl -E $SSL $SESSION/$ID?view=logs"
echo $cmd
$cmd

