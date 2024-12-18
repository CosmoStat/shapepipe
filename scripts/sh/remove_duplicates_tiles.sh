#!/usr/bin/env bash

count=0
for dir in tile_runs/*; do
  if compgen -G "$dir/output/run_sp_Mc_*" > /dev/null; then
    n=`ls -rdtl $dir/output/run_sp_Mc_* | wc -l`
    if [ "$n" != "2" ]; then
      ((n_remove=n-1))
      echo  $dir $n $n_remove
      rm -rf `ls -rdt1 $dir/output/run_sp_Mc_* | head -n $n_remove`
      (count=(count+1))
    fi
  fi
done

echo $count
