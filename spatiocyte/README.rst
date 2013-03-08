ecell4-spatiocyte
=========

$ cd ${ECELL3_PATH}
$ ./autogen.sh && ./configure --prefix=${PREFIX} --disable-python && make && make install
$ cd ${ECELL4_PATH}
$ PREFIX=~/local ECELL3_DM_PATH=${SPATIOCYTE_PATH} LIBRARY_PATH=${SPATIOCYTE_PATH} LD_LIBRARY_PATH=${SPATIOCYTE_PATH} CPLUS_INCLUDE_PATH=${PREFIX}/include/ecell-3.2:${SPATIOCYTE_PATH} ./install.sh spatiocyte
