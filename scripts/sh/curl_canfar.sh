#!/usr/bin/bash

# Usage
# ~/astro/repositories/github/shapepipe/scripts/sh/curl_canfar.sh 0.7 shapepipe/scripts/sh/init_canfar.sh 000.000

SSL=~/.ssl/cadcproxy.pem
NCORE=16
SESSION=https://ws-uv.canfar.net/skaha/v0/session
RESOURCES="ram=16&cores=$NCORE"
IMAGE=images.canfar.net/unions/shapepipe
NAME=shapepipe

# version of image on canfar, e.g. 0:7, 0:8
version=$1

# command on canfar, e.g. shapepipe/scripts/sh/init_run_canfar.sh
cmd=$2

# command line argument, e.g. 181.308
arg="$3 $NCORE"

#arg_an=`echo "$arg" | tr '_' 'X' | tr '.' 'X'`
arg_an=X

echo
echo "Start headless container"
echo "========================"
ID=`curl -E $SSL $SESSION?$RESOURCES -d "image=$IMAGE:$version" -d "name=${NAME}${arg_an}" -d "cmd=$cmd" --data-urlencode "args=$arg"`
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

