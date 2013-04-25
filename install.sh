#!/usr/bin/env bash

# PREFIX=/usr/local
# PREFIX=${HOME}/local
# PREFIX=
# SUBMODS=("bd" "gillespie")
SUBMODS=("bd" "gillespie" "ode" "egfrd" "spatiocyte")

CXXFLAGS="-g -Wall -Werror -Wno-uninitialized -O0 -DDEBUG" # enable debug mode

install_core()
{
    # install ecell4-core
    CXXFLAGS=${CXXFLAGS} ./waf distclean update --files="boost,doxygen" \
        configure --prefix=${PREFIX} build install
    return $?
}

install_submodule()
{
    # install a submodule
    if [ $# != 1 ]; then
        return 1
    fi
    cd $1
    LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${PREFIX}/lib \
        LIBRARY_PATH=${LIBRARY_PATH}:${PREFIX}/lib \
        CPLUS_INCLUDE_PATH=${CPLUS_INCLUDE_PATH}:${PREFIX}/include \
        CXXFLAGS=${CXXFLAGS} \
        ../waf distclean configure --prefix=${PREFIX} build install
    VAL=$?
    cd ..
    return ${VAL}
}

if [ "$PREFIX" == "" ]; then
    echo "\${PREFIX} is undefined."
    exit 1
fi

if [ $# == 0 ]; then
    TMP=("core" ${SUBMODS[@]})
else
    TMP=$@
fi

for SUBMOD in ${TMP[@]}
do
    if [ $SUBMOD == "core" ]; then
        install_core
    else
        install_submodule $SUBMOD
    fi
    if [ $? != 0 ]; then
        exit 1
    fi
done
