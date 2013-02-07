ecell4-egfrd
=========

$ cd ${EPDP_PATH}
$ ./autogen.sh && ./configure && make -j
$ ln -s _gfrd.so libgfrd.so
$ cd ${ECELL4_PATH}
$ LIBRARY_PATH=${EPDP_PATH} LD_LIBRARY_PATH=${EPDP_PATH} CPLUS_INCLUDE_PATH=${EPDP_PATH}:${ODEINT_PATH} ./install.sh egfrd
