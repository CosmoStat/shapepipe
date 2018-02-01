# !/bin/bash

# Set path to package module
export PACKAGE_DIR="$HOME/ShapePipe/modules/cfis_package"

# Run package
python ${PACKAGE_DIR}/cfis/cfis_SMP.py -d ${PACKAGE_DIR}/config -c package_config_smp.cfg
