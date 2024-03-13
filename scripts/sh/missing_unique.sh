#!/usr/bin/env bash

FILES=("summary/missing_job_128_ngmix_runner_*.txt")
temp="temp_temp.tmp"
temp2="temp_temp2.tmp"
out="missing_job_128_ngmix_runner_cut.txt"

i=0
for file in ${FILES[@]}; do

  echo $file $i

  if [ "$i" == "0" ]; then
    cp $file $temp
  else
    comm -12  <(sort $file) <(sort $temp) > $temp2
    cp $temp2 $temp
  fi

  wc $file $temp

  ((i=i+1))

done

mv $temp $out
rm $temp2
