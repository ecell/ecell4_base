#!/bin/bash -x

PYTHON_MAJOR_VERSION=$1
VTK_INCLUDE_PATH=/usr/include/vtk-5.8
WITH_VTK=0
WITH_HDF5=0

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
# rm -rf ${PREFIX}/lib/libecell4-*.so ${PREFIX}/include/ecell4;

set -e

cmake -DCMAKE_INSTALL_PREFIX=${PREFIX} -DWITH_HDF5=${WITH_HDF5} -DWITH_VTK=${WITH_VTK} .
make
# cmake -DCMAKE_INSTALL_PREFIX=${PREFIX} -DECELL4_ENABLE_PROFILING=1 .
# make VERBOSE=1
make test
make install

cd python

if [ "$PYTHON_MAJOR_VERSION" = "py2" ]; then
    # rm -rf build lib/ecell4/*.cpp
    mkdir -p ${PREFIX}/lib/python2.7/site-packages
    if [ "$(uname)" == "Darwin" ]; then
        LD_LIBRARY_PATH=${PREFIX}/lib PYTHONPATH=${PREFIX}/lib/python2.7/site-packages:/usr/local/lib/python2.7/dist-packages:${PYTHONPATH} python setup.py build_ext -L${PREFIX}/lib -I${PREFIX}/include:${VTK_INCLUDE_PATH} install --prefix=${PREFIX}
        PYTHONPATH=${PREFIX}/lib/python2.7/site-packages:/usr/local/lib/python2.7/dist-packages:${PYTHONPATH} LD_LIBRARY_PATH=${PREFIX}/lib python setup.py test
    else
        LD_LIBRARY_PATH=${PREFIX}/lib PYTHONPATH=${PREFIX}/lib/python2.7/site-packages:/usr/local/lib/python2.7/dist-packages:${PYTHONPATH} python2 setup.py build_ext -L${PREFIX}/lib -I${PREFIX}/include:${VTK_INCLUDE_PATH} install --prefix=${PREFIX}
        PYTHONPATH=${PREFIX}/lib/python2.7/site-packages:/usr/local/lib/python2.7/dist-packages:${PYTHONPATH} LD_LIBRARY_PATH=${PREFIX}/lib python2 setup.py test
    fi
elif [ "$PYTHON_MAJOR_VERSION" = "py3" ]; then
    # rm -rf build lib/ecell4/*.cpp
    mkdir -p ${PREFIX}/lib/python3.4/site-packages
    LD_LIBRARY_PATH=${PREFIX}/lib PYTHONPATH=${PREFIX}/lib/python3.4/site-packages:/usr/local/lib/python3.4/dist-packages:${PYTHONPATH} python3 setup.py build_ext -L${PREFIX}/lib -I${PREFIX}/include:${VTK_INCLUDE_PATH} install --prefix=${PREFIX}
    PYTHONPATH=${PREFIX}/lib/python3.4/site-packages:/usr/local/lib/python3.4/dist-packages:${PYTHONPATH} LD_LIBRARY_PATH=${PREFIX}/lib python3 setup.py test
fi
