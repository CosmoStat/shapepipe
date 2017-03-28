#!/bin/csh

# Set path to PPE modules
setenv PPE_DIR $HOME/pipeline/ShapePipe/modules/ppe_package
setenv PPE_CONFIG_DIR ${PPE_DIR}/config/

# Run PPE
python ${PPE_DIR}/ppe/ppe_SMP.py -d ${PPE_CONFIG_DIR} -c ppe_config_smp.cfg
