ecell4-egfrd
=========

$ cd ${EPDP_PATH}
$ ./autogen.sh && ./configure && make -j
$ cd ${ECELL4_PATH}
$ LIBRARY_PATH=${EPDP_PATH} LD_LIBRARY_PATH=${EPDP_PATH} CPLUS_INCLUDE_PATH=${EPDP_PATH} ./install.sh egfrd
