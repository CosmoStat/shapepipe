#!/bin/bash

inc_path="summary/special_job_128_ngmix_runner_6.txt"

cat $inc_path | perl -F[-.] -ane 'print "$F[1].$F[2]\n"' > inc.txt

for ID in `cat inc.txt`; do
  src=" tile_runs/$ID/output//run_sp_tile_ngmix_Ng1u"
  dst="tile_runs/$ID/output//run_sp_tile_ngmix_Ng1u_prev"

  if [ ! -e $dst ]; then

    # Previous run does not exist: move current to previous run
    mv $src $dst

  elif [ ! -e $src ]; then

    # Current does not exit (any more), skipping
    echo "$ID no current run found, skipping..."

  else

    size_src=$(stat -c %s $src/ngmix_runner/output/ngmix-*.fits)
    size_dst=$(stat -c %s $dst/ngmix_runner/output/ngmix-*.fits)

    echo $size_src $size_dst
    if ((size_src > size_dst)); then
      # Current output is larger than previous: overwrite
      rm -r $dst/ngmix_runner $dst/logs $dst/tmp
      mv $src/ngmix_runner $dst
      mv $src/logs $dst
      mv $src/tmp $dst
    else
      # Current output is not larger: remove current
      rm -r $src
    fi
  fi

done
