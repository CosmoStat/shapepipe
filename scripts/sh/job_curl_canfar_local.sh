# Usage:
# Edit the file "all.txt"
# screen
# bash run_curl.sh kind job

rm -f session_IDs.txt session_image_IDs.txt

script_local=curl_canfar_local.sh
version="1.0"
cmd_remote="shapepipe/scripts/sh/init_run_exclusive_canfar.sh"
N_SMP=1
kind="$1"
job="$2"

echo $kind $job
cat all.txt | xargs -n 1 -P 1 $script_local -v $version -c $cmd_remote -N $N_SMP -k $kind -j $job -e

