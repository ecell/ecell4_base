find_package(PythonInterp ${PY_VERSION})

set(PYTHON_VERSION ${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR})

if (PYTHONINTERP_FOUND)
    # set(ENV{LD_LIBRARY_PATH} ${BINARY_DIRECTORIES}:$ENV{LD_LIBRARY_PATH}})
    # set(ENV{PYTHONPATH} ${CMAKE_INSTALL_PREFIX}/lib/python${PYTHON_VERSION}/site-packages)

    set(SETUP_PY_IN ${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in)
    set(SETUP_PY    ${CMAKE_CURRENT_SOURCE_DIR}/setup.py)
    set(INIT_PY_IN  ${CMAKE_CURRENT_SOURCE_DIR}/lib/ecell4/__init__.py.in)
    set(INIT_PY     ${CMAKE_CURRENT_SOURCE_DIR}/lib/ecell4/__init__.py)
    set(OUTPUT      ${CMAKE_CURRENT_BINARY_DIR}/build/timestanp)

    configure_file(${SETUP_PY_IN} ${SETUP_PY})
    configure_file(${INIT_PY_IN} ${INIT_PY})

    if (NOT NO_SHARED)
        set(OPTS --prefer-shared)
    endif()

    add_custom_command(
        OUTPUT ${OUTPUT}
        COMMAND ${PYTHON_EXECUTABLE} ${SETUP_PY} build_ext -L${BINARY_DIRECTORIES} ${OPTS}
        COMMAND ${CMAKE_COMMAND} -E touch ${OUTPUT})
    add_custom_target(python ALL DEPENDS ${OUTPUT})
    add_test(
        NAME python_test
        COMMAND ${PYTHON_EXECUTABLE} ${SETUP_PY} test ${OPTS})
    install(CODE "
        set(ENV{PYTHONPATH} ${CMAKE_INSTALL_PREFIX}/lib/python${PYTHON_VERSION}/site-packages)
        set(ENV{LIBRARY_PATH} ${CMAKE_INSTALL_PREFIX}/lib)
        execute_process(COMMAND ${PYTHON_EXECUTABLE} ${SETUP_PY} install --prefix=${CMAKE_INSTALL_PREFIX} ${OPTS})")
endif()
