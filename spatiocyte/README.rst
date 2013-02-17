ecell4-spatiocyte
=========

Rename SpatiocyteStepper.so as libecell3-spatiocyte-stepper.so.

$ cd ${ECELL4_PATH}
$ PREFIX=~/local LD_LIBRARY_PATH=${SPATIOCYTE_PATH} CPLUS_INCLUDE_PATH=${PREFIX}/include/ecell-3.2:${PREFIX}/include/ecell-3.2/libecs:${SPATIOCYTE_PATH} ./install.sh spatiocyte
