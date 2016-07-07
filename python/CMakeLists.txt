find_package(PythonInterp ${PY_VERSION})

set(PYTHON_VERSION ${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR})

if (PYTHONINTERP_FOUND)
    set(ENV{PYTHONPATH} ${CMAKE_INSTALL_PREFIX}/lib/python${PYTHON_VERSION}/site-packages)
    set(ENV{INCLUDE_PATH})

    set(SETUP_PY_IN ${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in)
    set(SETUP_PY    ${CMAKE_CURRENT_BINARY_DIR}/setup.py)
    set(OUTPUT      ${CMAKE_CURRENT_BINARY_DIR}/build/timestanp)

    configure_file(${SETUP_PY_IN} ${SETUP_PY})

    add_custom_command(
        OUTPUT ${OUTPUT}
        COMMAND ${PYTHON_EXECUTABLE} ${SETUP_PY} build_ext
        COMMAND ${CMAKE_COMMAND} -E touch ${OUTPUT})
    add_custom_target(target ALL DEPENDS ${OUTPUT})
    add_test(
        NAME python_test
        COMMAND ${PYTHON_EXECUTABLE} ${SETUP_PY} test)
    install(
        CODE "execute_process(
        COMMAND ${PYTHON_EXECUTABLE} ${SETUP_PY} install --prefix=${CMAKE_INSTALL_PREFIX})")
endif()
