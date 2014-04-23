#!/usr/bin/env bash

# PREFIX=/usr/local
# PREFIX=${HOME}/local
# PREFIX=
# SUBMODS=("bd" "gillespie" "ode" "egfrd" "spatiocyte")
# PYTHONMODS=("core_python" "bd_python" "gillespie_python" "ode_python" "spatiocyte_python" "egfrd_python" "reaction_reader" "util_python")
SUBMODS=("bd" "gillespie" "ode" "lattice")
PYTHONMODS=("core_python" "bd_python" "gillespie_python" "ode_python")

CXXFLAGS="-g -Wall -Werror -Wno-uninitialized -O0 -DDEBUG" # enable debug mode
# WAFFLAGS="-v"

uninstall_all()
{
    rm -rf ${PREFIX}/include/ecell4
    rm -rf ${PREFIX}/lib/libecell4-*.so
    rm -rf ${PREFIX}/python2.7/site-packages/ecell4
}

install_core()
{
    # install ecell4-core
    cd core
    CXXFLAGS=${CXXFLAGS} python ../waf distclean update --files="boost,doxygen" \
        configure --prefix=${PREFIX} ${WAFFLAGS} install
    cd ..
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
        python ../waf distclean configure --prefix=${PREFIX} ${WAFFLAGS} install
    VAL=$?
    cd ..
    return ${VAL}
}

if [ "$PREFIX" == "" ]; then
    echo "\${PREFIX} is undefined."
    exit 1
fi

if [ $# == 0 ]; then
    TMP=("uninstall" "core" ${SUBMODS[@]} ${PYTHONMODS[@]})
else
    TMP=$@
fi

for SUBMOD in ${TMP[@]}
do
    if [ $SUBMOD == "uninstall" ]; then
        uninstall_all
    elif [ $SUBMOD == "core" ]; then
        install_core
    else
        install_submodule $SUBMOD
    fi
    if [ $? != 0 ]; then
        exit 1
    fi
done
