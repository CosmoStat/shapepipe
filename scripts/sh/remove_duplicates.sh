#!/usr/bin/env bash

count=0
for dir in exp_runs/*; do
  if compgen -G "$dir/output/run_sp_exp_Sx*" > /dev/null; then
    n=`ls -rdtl $dir/output/run_sp_exp_Sx* | wc -l`
    #echo $dir $n $count
    if [ "$n" == "2" ]; then
      n_remove=$((n-1))
      rm -rf `ls -rdt1 $dir/output/run_sp_exp_Sx* | head -n $n_remove`
      echo rm -rf `ls -rdt1 $dir/output/run_sp_exp_Sx* | head -n $n_remove`
      count=$((count+1))
    fi
  fi
done

echo $count
