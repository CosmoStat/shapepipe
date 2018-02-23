#! /bin/bash

# Set path to package module
export PACKAGE_DIR="/Users/aguinot/Documents/pipeline/pipeline_dev/MakeVignet_package"

# Run package
python ${PACKAGE_DIR}/MakeVignet/MakeVignet_SMP.py -d ${PACKAGE_DIR}/config -c package_config_smp.cfg
