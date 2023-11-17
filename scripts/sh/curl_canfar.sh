#!/usr/bin/bash

# Usage
# cat tile_num.txt | xargs -n 1 -P 1 ~/astro/repositories/github/shapepipe/scripts/sh/curl_canfar.sh 0.7 shapepipe/scripts/sh/init_canfar.sh

SSL=~/.ssl/cadcproxy.pem
SESSION=https://ws-uv.canfar.net/skaha/v0/session
IMAGE=images.canfar.net/unions/shapepipe
NAME=shapepipe

version=$1
cmd=$2
arg=$3

arg_an=`echo "$arg" | tr '_' 'X' | tr '.' 'X'`
arg_an=X

echo "Start headless container"
echo "========================"
ID=`curl -E $SSL $SESSION -d "image=$IMAGE:$version" -d "name=${NAME}${arg_an}" -d "cmd=$cmd" --data-urlencode "args=$arg"`
#curl -E $SSL $SESSION -d \"name=$NAME\" -d \"image=images.canfar.net/unions/shapepipe\:0.7\" -d \"cmd=$2\" --data-urlencode \"args=$3\"
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

