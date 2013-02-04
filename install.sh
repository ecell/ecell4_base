#!/usr/bin/env bash

# PREFIX=/usr/local
# PREFIX=${HOME}/local
PREFIX=
SUBMODS=("bd" "gillespie")

if [ "$PREFIX" == "" ]; then
    echo "\${PREFIX} is undefined."
    exit 1
fi

# install ecell4-core
./waf distclean update --files="boost,doxygen" configure --prefix=${PREFIX} build install
if [ $? != 0 ]; then
    exit 1
fi

# install submodules
for SUBMOD in ${SUBMODS[@]}
do
    cd $SUBMOD
    LD_LIBRARY_PATH=${PREFIX}/lib LIBRARY_PATH=${PREFIX}/lib \
        CPLUS_INCLUDE_PATH=${PREFIX}/include \
        ../waf distclean configure --prefix=${PREFIX} build install
    cd ..
    if [ $? != 0 ]; then
        exit 1
  fi
done
