find_package(PythonInterp ${PY_VERSION})

if (PYTHONINTERP_FOUND)
    if (DEFINED PYTHON_VERSION_MAJOR AND DEFINED PYTHON_VERSION_MINOR)
        set(PYTHON_VERSION ${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR})
    else (DEFINED PYTHON_VERSION_MAJOR AND DEFINED PYTHON_VERSION_MINOR)
        execute_process(COMMAND ${PYTHON_EXECUTABLE} --version
            OUTPUT_VARIABLE PYTHON_VERSION_RAW
            ERROR_VARIABLE PYTHON_VERSION_RAW)
        string(REPLACE "\n" "" PYTHON_VERSION_RAW "${PYTHON_VERSION_RAW}")
        string(REGEX MATCH "[0-9]\\.[0-9]" PYTHON_VERSION "${PYTHON_VERSION_RAW}")
    endif (DEFINED PYTHON_VERSION_MAJOR AND DEFINED PYTHON_VERSION_MINOR)

    # set(ENV{LD_LIBRARY_PATH} ${BINARY_DIRECTORIES}:$ENV{LD_LIBRARY_PATH}})
    # set(ENV{PYTHONPATH} ${CMAKE_INSTALL_PREFIX}/lib/python${PYTHON_VERSION}/site-packages)

    get_directory_property(INCLUDE_PATH INCLUDE_DIRECTORIES)
    set(_TMP_RESULT "")
    set(_GLUE "") # effective glue is empty at the beginning
    foreach(SINGLE_PATH ${INCLUDE_PATH})
      set(_TMP_RESULT "${_TMP_RESULT}${_GLUE}${SINGLE_PATH}")
      set(_GLUE ":") # delimiter
    endforeach()
    set(INCLUDE_PATH ${_TMP_RESULT})
    # message("${INCLUDE_PATH}")
    message("XXX: CMAKE_CURRENT_SOURCE_DIR is '${CMAKE_CURRENT_SOURCE_DIR}'.")

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
        COMMAND ${PYTHON_EXECUTABLE} ${SETUP_PY} bdist_wheel
        COMMAND ${CMAKE_COMMAND} -E touch ${OUTPUT})
    add_custom_target(python ALL DEPENDS ${OUTPUT})
    add_test(
        NAME python_test
        COMMAND ${PYTHON_EXECUTABLE} ${SETUP_PY} test ${OPTS})
    install(DIRECTORY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/python${PYTHON_VERSION}/site-packages)
    install(CODE "
        set(ENV{PYTHONPATH} ${CMAKE_INSTALL_PREFIX}/lib/python${PYTHON_VERSION}/site-packages)
        set(ENV{LIBRARY_PATH} ${CMAKE_INSTALL_PREFIX}/lib)
        execute_process(COMMAND ${PYTHON_EXECUTABLE} ${SETUP_PY} install --prefix=${CMAKE_INSTALL_PREFIX} ${OPTS})")
endif()
