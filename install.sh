#!/bin/bash -x

CURDIR=$(cd $(dirname $0); pwd)
# PREFIX=/usr/local
# PREFIX=${HOME}/local
# PREFIX=${CURDIR}/local
# PREFIX=

# make clean; rm -rf ${PREFIX}; rm CMakeCache.txt
# rm ecell4/egfrd/SphericalBesselTable.hpp ecell4/egfrd/CylindricalBesselTable.hpp

set -e

cd ecell4/egfrd/tablegen
cmake .
make
cp SphericalBesselTable.hpp CylindricalBesselTable.hpp ..
cd ../../..

cmake -DCMAKE_INSTALL_PREFIX=${PREFIX} .
make
# cmake -DCMAKE_INSTALL_PREFIX=${PREFIX} -DECELL4_ENABLE_PROFILING=1 .
# make VERBOSE=1
make test
make install

cd python

# rm -rf build lib/ecell4/*.cpp
mkdir -p ${PREFIX}/lib/python2.7/site-packages
LD_LIBRARY_PATH=${PREFIX}/lib PYTHONPATH=${PREFIX}/lib/python2.7/site-packages:/usr/local/lib/python2.7/dist-packages:${PYTHONPATH} python setup.py build_ext -L${PREFIX}/lib -I${PREFIX}/include install --prefix=${PREFIX}
PYTHONPATH=${PREFIX}/lib/python2.7/site-packages:/usr/local/lib/python2.7/dist-packages:${PYTHONPATH} LD_LIBRARY_PATH=${PREFIX}/lib python setup.py test

# rm -rf build lib/ecell4/*.cpp
# mkdir -p ${PREFIX}/lib/python3.4/site-packages
# LD_LIBRARY_PATH=${PREFIX}/lib PYTHONPATH=${PREFIX}/lib/python3.4/site-packages:/usr/local/lib/python3.4/dist-packages:${PYTHONPATH} python3 setup.py build_ext -L${PREFIX}/lib -I${PREFIX}/include install --prefix=${PREFIX}
# PYTHONPATH=${PREFIX}/lib/python3.4/site-packages:/usr/local/lib/python3.4/dist-packages:${PYTHONPATH} LD_LIBRARY_PATH=${PREFIX}/lib python3 setup.py test
