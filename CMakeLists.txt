cmake_minimum_required(VERSION 2.6)
# cmake_minimum_required(VERSION 3.0)
project(ecell4)
set(PROJECT_VERSION 1.0.dev2)

# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11 -fext-numeric-literals")  #XXX: c++11

if(ECELL4_ENABLE_PROFILING)
    add_definitions(-pg -O0)
    set(CMAKE_EXE_LINKER_FLAGS "-pg -O0 ${CMAKE_EXE_LINKER_FLAGS}")
    set(CMAKE_SHARED_LINKER_FLAGS "-pg -O0 ${CMAKE_SHARED_LINKER_FLAGS}")
    set(CMAKE_MODULE_LINKER_FLAGS "-pg -O0 ${CMAKE_MODULE_LINKER_FLAGS}")
else(ECELL4_ENABLE_PROFILING)
    add_definitions(-O2 -Wall)
    set(CMAKE_EXE_LINKER_FLAGS "-O2 -Wall ${CMAKE_EXE_LINKER_FLAGS}")
    set(CMAKE_SHARED_LINKER_FLAGS "-O2 -Wall ${CMAKE_SHARED_LINKER_FLAGS}")
    set(CMAKE_MODULE_LINKER_FLAGS "-O2 -Wall ${CMAKE_MODULE_LINKER_FLAGS}")
    #XXX:STRICT: add_definitions(-O2 -Wall -Werror)
    #XXX:STRICT: set(CMAKE_EXE_LINKER_FLAGS "-O2 -Wall -Werror ${CMAKE_EXE_LINKER_FLAGS}")
    #XXX:STRICT: set(CMAKE_SHARED_LINKER_FLAGS "-O2 -Wall -Werror ${CMAKE_SHARED_LINKER_FLAGS}")
    #XXX:STRICT: set(CMAKE_MODULE_LINKER_FLAGS "-O2 -Wall -Werror ${CMAKE_MODULE_LINKER_FLAGS}")
endif(ECELL4_ENABLE_PROFILING)
include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_BINARY_DIR})
enable_testing()

find_package(VTK QUIET)
if(VTK_FOUND)
  include(${VTK_USE_FILE})
  set(WITH_VTK 1)
  set(HAVE_VTK 1)
else()
  set(WITH_VTK 0)
endif()

message("Looking for HDF5")
find_package(HDF5 COMPONENTS C CXX HL REQUIRED)
# find_package(HDF5 COMPONENTS C CXX HL QUIET)
if (HDF5_FOUND)
  message("HDF5 was found ... ${HDF5_VERSION}")
  include_directories(${HDF5_INCLUDE_DIRS})
  set(WITH_HDF5 1)
else()
  message("Could NOT find HDF5")
  set(HDF5_LIBRARIES)
  set(WITH_HDF5 0)
endif()

find_package(Boost COMPONENTS regex unit_test_framework)
if(NOT DEFINED Boost_VERSION)
    message(FATAL_ERROR
        "ecell4 requires Boost C++ Library. "
        "If you have already installed it, try `-DBOOST_ROOT=/path/to/boost`")
else()
    include_directories(${Boost_INCLUDE_DIRS})
    link_directories(${Boost_LIBRARY_DIRS})
endif()

if(Boost_REGEX_FOUND)
  set(HAVE_BOOST_REGEX 1)
else()
  set(HAVE_BOOST_REGEX 0)
endif()

find_library(GSL_LIBRARIES gsl)
find_library(GSL_CBLAS_LIBRARIES gslcblas)

include(CheckIncludeFileCXX)
# set(CMAKE_REQUIRED_INCLUDES "/usr/include")
include(CheckCXXSourceCompiles)

check_include_file_cxx(unordered_map HAVE_UNORDERED_MAP)
check_cxx_source_compiles("#include <tr1/unordered_map>
int main() { std::tr1::unordered_map<int, int> a, b(a); return 0; }"
    HAVE_TR1_UNORDERED_MAP)
check_include_file_cxx("boost/unordered_map.hpp" HAVE_BOOST_UNORDERED_MAP_HPP)

check_cxx_source_compiles("#include <tr1/functional>
int main() { std::tr1::hash<int>(); return 0; }"
    HAVE_TR1_FUNCTIONAL)
# check_include_file_cxx(functional HAVE_STD_HASH)
check_cxx_source_compiles("#include <functional>
int main() { std::hash<int>(); return 0; }"
    HAVE_STD_HASH)
check_include_file_cxx("boost/functional/hash.hpp" HAVE_BOOST_FUNCTIONAL_HASH_HPP)

check_include_file_cxx(chrono HAVE_CHRONO)

include_directories(${CMAKE_BINARY_DIR})

# if (NOT NO_SHARED)
  add_subdirectory(greens_functions)  #XXX: submodule
# endif()

find_package(PythonInterp ${PY_VERSION})

execute_process(
    COMMAND ${PYTHON_EXECUTABLE} -c "from distutils import sysconfig; print(sysconfig.get_python_inc())"
    OUTPUT_VARIABLE PYTHON_INCLUDE_DIR)

include_directories(${PYTHON_INCLUDE_DIR})
unset(PYTHON_INCLUDE_DIR)

add_subdirectory(ecell4)
add_subdirectory(python)
