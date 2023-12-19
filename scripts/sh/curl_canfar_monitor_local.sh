#!/usr/bin/bash

# --false-start --tcp-fastopen  faster

# --http2   slower
# --tlsv1.2 slower?
# --insecure ?

# -H "Accept-Encoding: gzip"    faster?

SSL=~/.ssl/cadcproxy.pem
SESSION=https://ws-uv.canfar.net/skaha/v0/session

type=$1

echo "type=$type"

while [ 1 ]; do
  session_ID=`tail -n 1 session_IDs.txt`
  cmd="curl -E $SSL $SESSION/$session_ID?view=$type"
  echo $cmd
  $cmd
done
