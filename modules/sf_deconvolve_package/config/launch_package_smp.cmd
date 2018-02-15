#!/bin/csh

# Set path to PPE modules
setenv PACKAGE_DIR $HOME/pipeline/ShapePipe/modules/sf_deconvolve_package
setenv PACKAGE_CONFIG_DIR ${PACKAGE_DIR}/config/

# Run PPE
python ${PACKAGE_DIR}/sf_deconvolve_package/sf_deconvolve_SMP.py -d ${PACKAGE_CONFIG_DIR} -c package_config_smp.cfg
