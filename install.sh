#!/bin/bash -x

PYTHON_MAJOR_VERSION=$1

if [ "${PREFIX-UNDEF}" = "UNDEF" ]; then
    if [ "$PREFIX" = "" ]; then
        PREFIX=/usr/local
        # PREFIX=${HOME}/local
        # CURDIR=$(cd $(dirname $0); pwd)
        # PREFIX=${CURDIR}/local
        # PREFIX=
    fi
fi

# make clean; rm CMakeCache.txt
# rm ecell4/egfrd/SphericalBesselTable.hpp ecell4/egfrd/CylindricalBesselTable.hpp
# rm -rf ${PREFIX}/lib/libecell4-*.so ${PREFIX}/include/ecell4;

set -e

if [ ! -f ecell4/egfrd/SphericalBesselTable.hpp -o ! -f ecell4/egfrd/CylindricalBesselTable.hpp ]; then
    cd ecell4/egfrd/tablegen
    cmake .
    make
    cp SphericalBesselTable.hpp CylindricalBesselTable.hpp ..
    cd ../../..
fi

cmake -DCMAKE_INSTALL_PREFIX=${PREFIX} .
make
# cmake -DCMAKE_INSTALL_PREFIX=${PREFIX} -DECELL4_ENABLE_PROFILING=1 .
# make VERBOSE=1
make test
make install

cd python

if [ "$PYTHON_MAJOR_VERSION" = "py2" ]; then
    # rm -rf build lib/ecell4/*.cpp
    mkdir -p ${PREFIX}/lib/python2.7/site-packages
    LD_LIBRARY_PATH=${PREFIX}/lib PYTHONPATH=${PREFIX}/lib/python2.7/site-packages:/usr/local/lib/python2.7/dist-packages:${PYTHONPATH} python setup.py build_ext -L${PREFIX}/lib -I${PREFIX}/include install --prefix=${PREFIX}
    PYTHONPATH=${PREFIX}/lib/python2.7/site-packages:/usr/local/lib/python2.7/dist-packages:${PYTHONPATH} LD_LIBRARY_PATH=${PREFIX}/lib python setup.py test
elif [ "$PYTHON_MAJOR_VERSION" = "py3" ]; then
    # rm -rf build lib/ecell4/*.cpp
    mkdir -p ${PREFIX}/lib/python3.4/site-packages
    LD_LIBRARY_PATH=${PREFIX}/lib PYTHONPATH=${PREFIX}/lib/python3.4/site-packages:/usr/local/lib/python3.4/dist-packages:${PYTHONPATH} python3 setup.py build_ext -L${PREFIX}/lib -I${PREFIX}/include install --prefix=${PREFIX}
    PYTHONPATH=${PREFIX}/lib/python3.4/site-packages:/usr/local/lib/python3.4/dist-packages:${PYTHONPATH} LD_LIBRARY_PATH=${PREFIX}/lib python3 setup.py test
fi
