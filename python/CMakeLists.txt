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

    set(ECELL4_SHARED_DIRS ${PROJECT_BINARY_DIR}/greens_functions:${ECELL4_SHARED_DIRS})
    include_directories(${PROJECT_SOURCE_DIR}/greens_functions ${PROJECT_BINARY_DIR}/greens_functions)

    get_directory_property(INCLUDE_PATH INCLUDE_DIRECTORIES)
    set(_TMP_RESULT "")
    set(_GLUE "") # effective glue is empty at the beginning
    foreach(SINGLE_PATH ${INCLUDE_PATH})
      set(_TMP_RESULT "${_TMP_RESULT}${_GLUE}${SINGLE_PATH}")
      set(_GLUE ";") # delimiter
    endforeach()
    set(INCLUDE_PATH ${_TMP_RESULT})
    message("XXX: INCLUDE_PATH '${INCLUDE_PATH}'")

    set(_TMP_RESULT "")
    set(_GLUE "") # effective glue is empty at the beginning
    foreach(SINGLE_PATH ${HDF5_CXX_LIBRARIES})
      get_filename_component(HDF5_LIBRARY_DIRS ${SINGLE_PATH} DIRECTORY)
      set(_TMP_RESULT "${_TMP_RESULT}${_GLUE}${HDF5_LIBRARY_DIRS}")
      set(_GLUE ";") # delimiter
    endforeach()
    set(HDF5_LIBRARY_DIRS ${_TMP_RESULT})

    # message("XXX: CMAKE_CURRENT_SOURCE_DIR '${CMAKE_CURRENT_SOURCE_DIR}'")
    # message("XXX: HDF5_LIBRARY_DIRS '${HDF5_LIBRARY_DIRS}'")
    # message("XXX: VTK_LIBRARIES '${VTK_LIBRARIES}'")
    # message("XXX: VTK_LIBRARY_DIRS '${VTK_LIBRARY_DIRS}'")
    # message("XXX: ECELL4_SHARED_DIRS '${ECELL4_SHARED_DIRS}'")

    set(SETUP_PY_IN ${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in)
    # set(SETUP_PY    ${CMAKE_CURRENT_SOURCE_DIR}/setup.py)
    set(SETUP_PY    ${CMAKE_CURRENT_BINARY_DIR}/setup.py)
    set(INIT_PY_IN  ${CMAKE_CURRENT_SOURCE_DIR}/lib/ecell4/__init__.py.in)
    set(INIT_PY     ${CMAKE_CURRENT_SOURCE_DIR}/lib/ecell4/__init__.py)
    set(OUTPUT      ${CMAKE_CURRENT_BINARY_DIR}/build/timestanp)

    configure_file(${SETUP_PY_IN} ${SETUP_PY})
    configure_file(${INIT_PY_IN} ${INIT_PY})

    if (NO_SHARED)
        add_custom_command(
            OUTPUT ${OUTPUT}
            COMMAND ${PYTHON_EXECUTABLE} ${SETUP_PY} build_ext
            COMMAND ${PYTHON_EXECUTABLE} ${SETUP_PY} bdist_wheel
            COMMAND ${CMAKE_COMMAND} -E touch ${OUTPUT})
    else()
        add_custom_command(
            OUTPUT ${OUTPUT}
            COMMAND ${PYTHON_EXECUTABLE} ${SETUP_PY} build_ext -L${ECELL4_SHARED_DIRS} --prefer-shared
            COMMAND ${CMAKE_COMMAND} -E touch ${OUTPUT})
    endif()

    add_custom_target(python ALL DEPENDS ${OUTPUT})
    add_test(
        NAME python_test
        COMMAND ${PYTHON_EXECUTABLE} ${SETUP_PY} test ${OPTS})
    install(DIRECTORY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/python${PYTHON_VERSION}/site-packages)

    if (NO_SHARED)
        install(CODE "
            set(ENV{PYTHONPATH} ${CMAKE_INSTALL_PREFIX}/lib/python${PYTHON_VERSION}/site-packages)
            set(ENV{LIBRARY_PATH} ${CMAKE_INSTALL_PREFIX}/lib)
            execute_process(COMMAND ${PYTHON_EXECUTABLE} ${SETUP_PY} install --prefix=${CMAKE_INSTALL_PREFIX})")
    else()
        install(CODE "
            set(ENV{PYTHONPATH} ${CMAKE_INSTALL_PREFIX}/lib/python${PYTHON_VERSION}/site-packages)
            set(ENV{LIBRARY_PATH} ${CMAKE_INSTALL_PREFIX}/lib)
            execute_process(COMMAND ${PYTHON_EXECUTABLE} ${SETUP_PY} install --prefix=${CMAKE_INSTALL_PREFIX} --prefer-shared)")
    endif()
    #     execute_process(COMMAND ${PYTHON_EXECUTABLE} ${SETUP_PY} install --prefix=${CMAKE_INSTALL_PREFIX} ${OPTS})")
endif()
