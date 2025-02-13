#!/usr/bin/env bash

# -H "Accept-Encoding: gzip"    faster?

SSL=~/.ssl/cadcproxy.pem
SESSION=https://ws-uv.canfar.net/skaha/v0/session

type=$1

echo "type=$type"

for session_ID in `cat ids.txt`; do
  cmd="curl -E $SSL $SESSION/$session_ID?view=$type"
  echo $cmd
  $cmd
done

exit 0

while [ 1 ]; do
  session_ID=`tail -n 1 ids.txt`
  cmd="curl -E $SSL $SESSION/$session_ID?view=$type"
  echo $cmd
  $cmd
done
