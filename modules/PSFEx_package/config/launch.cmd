#! /bin/bash

# Set path to package module
export PACKAGE_DIR="$HOME/stage_cea/gitlab/new/ShapePipe/modules/PSFEx_package"

# Run package
python ${PACKAGE_DIR}/PSFEx/PSFEx_SMP.py -d ${PACKAGE_DIR}/config -c package_config_smp.cfg
