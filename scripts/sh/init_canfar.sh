#!/bin/bash

echo "start init canfar"

<<<<<<< HEAD
#echo init_canfar > ~/init_canfar.log
#date >> ~/init_canfar.log
=======
echo init_canfar > ~/init_canfar.log
date >> ~/init_canfar.log
>>>>>>> origin/exclusive

. /opt/conda/etc/profile.d/conda.sh

conda activate shapepipe

echo "end init canfar"

