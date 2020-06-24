#!/usr/bin/env bash

# Name: canfar_rm_results.sh
# Description: Removes result FITS files on vos from canfar run
# Author: Martin Kilbinger <martin.kilbinger@cea.fr>
# Date: 06/2020
# Package: shapepipe


# The shapepipe python virtual environment needs to 
# be active to run this script.

TILE_ARR=($@)

list_path=ID_list.txt

for TILE in ${TILE_ARR[@]}; do
    vls vos:cfis/cosmostat/kilbinger/results/*$TILE* > $list_path
    for i in `cat $list_path`; do
        echo "vrm -v --quick vos:cfis/cosmostat/kilbinger/results/$i"
    done
done
