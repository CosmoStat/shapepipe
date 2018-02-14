# !/bin/bash

# Set path to package module
export PACKAGE_DIR="$HOME/ShapePipe/modules/shapelens_package"

# Run package
python ${PACKAGE_DIR}/shapelens/shapelens_SMP.py -d ${PACKAGE_DIR}/config -c package_config_smp.cfg
