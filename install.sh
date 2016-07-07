#!/bin/bash -x

usage_exit() {
    echo "Usage: $0 [--[no]python2] [--[no]python3] [--nopython] [--[no]test] [--[no]vtk] [--[no]hdf5] [--clean] [--prefix=PREFIX] [-h] [--help] [py2|py3]" 1>&2
    exit 1
}

WITH_PYTHON2=1
WITH_PYTHON3=0
WITH_TEST=1
VTK_INCLUDE_PATH=/usr/include/vtk-5.8
WITH_VTK=0
WITH_HDF5=0
CLEANUP=0
NO_BESSEL_TABLE=0

if [ "${PREFIX-UNDEF}" = "UNDEF" ]; then
    if [ "$PREFIX" = "" ]; then
        PREFIX=/usr/local
        # PREFIX=${HOME}/local
        # CURDIR=$(cd $(dirname $0); pwd)
        # PREFIX=${CURDIR}/local
        # PREFIX=
    fi
fi

while getopts :h-: OPT
do
    case $OPT in
        -)
            case "$OPTARG" in
                python2) WITH_PYTHON2=1;;
                python3) WITH_PYTHON3=1;;
                nopython2) WITH_PYTHON2=0;;
                nopython3) WITH_PYTHON3=0;;
                nopython)
                    WITH_PYTHON2=0;
                    WITH_PYTHON3=0;
                    ;;
                test) WITH_TEST=1;;
                notest) WITH_TEST=0;;
                vtk) WITH_VTK=1;;
                novtk) WITH_VTK=0;;
                hdf5) WITH_HDF5=1;;
                nohdf5) WITH_HDF5=0;;
                prefix=*) PREFIX=${OPTARG#*=};;
                nobessel) NO_BESSEL_TABLE=1;;
                clean) CLEANUP=1;;
                help) usage_exit;;
                *) usage_exit;;
            esac;;
        h) usage_exit;;
        \?) usage_exit;;
    esac
done

shift $((OPTIND - 1))

PREFIX=$(cd $(dirname $PREFIX) && pwd)/$(basename $PREFIX)

PYTHON_MAJOR_VERSION=$1
if [ "$PYTHON_MAJOR_VERSION" = "py2" ]; then
    WITH_PYTHON2=1
    WITH_PYTHON3=0
elif [ "$PYTHON_MAJOR_VERSION" = "py3" ]; then
    WITH_PYTHON2=0
    WITH_PYTHON3=1
else
    echo "No Python major version was specified: WITH_PYTHON2=${WITH_PYTHON2} WITH_PYTHON3=${WITH_PYTHON3}"
fi

if [ $CLEANUP != 0 ]; then
    make clean; rm CMakeCache.txt
    rm -rf ${PREFIX}/lib/libecell4-*.so ${PREFIX}/include/ecell4;
fi

set -e

cmake \
  -DCMAKE_INSTALL_PREFIX=${PREFIX} \
  -DWITH_HDF5=${WITH_HDF5} \
  -DWITH_VTK=${WITH_VTK} \
  -DNO_BESSEL_TABLE=${NO_BESSEL_TABLE} .
make
# cmake -DCMAKE_INSTALL_PREFIX=${PREFIX} -DECELL4_ENABLE_PROFILING=1 .
# make VERBOSE=1
if [ $WITH_TEST != 0 ]; then
    make test
fi
make install

cd python

PYTHON=python
INCLUDE_PATH=${PREFIX}/include:${VTK_INCLUDE_PATH}

if [ $WITH_PYTHON2 != 0 ]; then
    if [ $CLEANUP != 0 ]; then
        rm -rf build lib/ecell4/*.cpp
    fi

    if [ "$(uname)" != "Darwin" ]; then
        PYTHON=python2
        if [ $WITH_HDF5 != 0 ]; then
          INCLUDE_PATH=/usr/include/hdf5/serial/:${INCLUDE_PATH}
        fi
        BUILD_OPT="--prefer-shared"
    fi
elif [ $WITH_PYTHON3 != 0 ]; then
    if [ $CLEANUP != 0 ]; then
        rm -rf build lib/ecell4/*.cpp
    fi

    PYTHON=python3
    BUILD_OPT="--prefer-shared"
else
    exit 0
fi

python_ver() {
    $1 --version | awk '{print $2;}' | cut -d'.' -f 1,2
}

PYTHON_VERSION=$(python_ver ${PYTHON})
DEST_DIR=${PREFIX}/lib/python${PYTHON_VERSION}/site-packages
mkdir -p ${DEST_DIR}

export PYTHONPATH=${DEST_DIR}:/usr/local/lib/python${PYTHON_VERSION}/dist-packages:${PYTHONPATH}
export LD_LIBRARY_PATH=${PREFIX}/lib:${LD_LIBRARY_PATH}

${PYTHON} setup.py build_ext -L${PREFIX}/lib -I${INCLUDE_PATH} ${BUILD_OPT}
${PYTHON} setup.py install --prefix=${PREFIX}
${PYTHON} setup.py test
