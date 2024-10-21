#!/usr/bin/env bash

# Name: stats_jobs_canfar.sh
# Author: Martin Kilbinger <martin.kilbinger@cea.fr>
# Description: Handles headless jobs on canfar


# Global variables

## Temporary files
tmpfile_jobs="jobinfo.txt"
tmpfile_ids="ids.txt"
tmpfile_running="jobs_running.txt"

## curl options
SSL=~/.ssl/cadcproxy.pem
SESSION=https://ws-uv.canfar.net/skaha/v0/session


# Command line arguments

## Default values
mode="count"

## Help string
usage="Usage: $(basename "$0") [OPTIONS]
\n\nOptions:\n
   -h\tthis message\n
   -m, --mode MODE\n
   \tmode, allowed are 'count' (default), 'delete'\n
"

## Parse command line
while [ $# -gt 0 ]; do
  case "$1" in
    -h)
      echo -ne $usage
      exit 0
      ;;
    -m|--mode)
      mode="$2"
      shift
      ;;
  esac
  shift
done

## Check options
case $mode in
  "count"|"delete")
    # valid option
    ;;
  *)
    echo "Invalid mode $mode"
    exit 1
    ;;
esac


# Main program

# Get all instances
#echo curl -E $SSL $SESSION
curl -E $SSL $SESSION &> /dev/null > $tmpfile_jobs
res=$?

if [ "$res" == "0" ]; then
  # Number of jobs
  n_headless=`cat $tmpfile_ids | grep Running | wc -l`
else
    # Failure: set to very high number
    n_headless=10000
fi

# Get headless job IDs
#cat $tmpfile_jobs | grep headless -B 4 -A 12 | grep \"id | perl -F\" -ane 'print "$F[3]\n"' > $tmpfile_ids
cat $tmpfile_jobs | grep headless -B 4 -A 2 | grep Running -A 1 > $tmpfile_ids


# Get running job info
cat $tmpfile_ids | grep name | perl -F\- -ane 'chomp; $F[4] =~ s/[",]//g; print "$F[3].$F[4]"' > $tmpfile_running


if [ "$mode" == "count" ]; then

  echo $n_headless

elif [ "$mode" == "delete" ]; then

  echo -n "Delete $n_headless jobs? [y|n] "
  read answer
  if [ "$answer" == "y" ]; then
    cat $tmpfile_jobs | grep headless -B 34 -A 6 | grep Running -A 34 | grep id | grep -v user | perl -F\"  -ane 'print "$F[3]\n"' > $tmpfile_ids
    for ID in `cat $tmpfile_ids`; do
      echo $ID
      # Delete headless jobs
      echo "curl -X DELETE -E $SSL $SESSION/$ID"
      curl -X DELETE -E $SSL $SESSION/$ID
      echo $?
    done
  fi

fi


# Remove temporary files
#rm -f $tmpfile_jobs $tmpfile_ids
