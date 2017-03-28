#!/bin/csh

# Set path to PPE modules
setenv PPE_DIR $HOME/pipeline/ShapePipe/modules/ppe_package

# Run PPE
python ${PPE_DIR}/ppe/ppe_SMP.py -d ${PPE_DIR}/config/ -c ppe_config_smp.cfg