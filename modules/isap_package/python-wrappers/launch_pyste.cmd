#!/bin/csh
cd isap-great3-wrappers/tools
pyste -I${EUCLID_DIR}/code/cea-epfl/branches/isap_package/src/libtools -I${EUCLID_DIR}/code/cea-epfl/branches/isap_package/src/libsparse1d -I${EUCLID_DIR}/code/cea-epfl/branches/isap_package/src/libsparse2d --module=tools tools.pyste
