# Usage:
# Edit the file "all.txt"
# screen
# bash run_curl.sh kind job

rm -f session_IDs.txt session_image_IDs.txt

script_local=$HOME/astro/repositories/github/shapepipe/scripts/sh/curl_canfar_local.sh
version="0.9"
cmd_remote="shapepipe/scripts/sh/init_run_exclusive_canfar.sh"
kind="$1"
job="$2"

echo $kind $job
cat all.txt | xargs -n 1 -P 10 $script_local -v $version -c $cmd_remote -k $kind -j $job -e

