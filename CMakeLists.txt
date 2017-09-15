cmake_minimum_required(VERSION 2.6)
# cmake_minimum_required(VERSION 3.0)
project(ecell4)
set(PROJECT_VERSION 4.1.2)

if(ECELL4_ENABLE_PROFILING)
    add_definitions(-pg -O0)
    set(CMAKE_EXE_LINKER_FLAGS "-pg -O0 ${CMAKE_EXE_LINKER_FLAGS}")
    set(CMAKE_SHARED_LINKER_FLAGS "-pg -O0 ${CMAKE_SHARED_LINKER_FLAGS}")
    set(CMAKE_MODULE_LINKER_FLAGS "-pg -O0 ${CMAKE_MODULE_LINKER_FLAGS}")
else(ECELL4_ENABLE_PROFILING)
    add_definitions(-O2)
    set(CMAKE_EXE_LINKER_FLAGS "-O2 ${CMAKE_EXE_LINKER_FLAGS}")
    set(CMAKE_SHARED_LINKER_FLAGS "-O2 ${CMAKE_SHARED_LINKER_FLAGS}")
    set(CMAKE_MODULE_LINKER_FLAGS "-O2 ${CMAKE_MODULE_LINKER_FLAGS}")
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

find_package(HDF5 COMPONENTS C CXX HL QUIET)
if (HDF5_FOUND)
  include_directories(${HDF5_INCLUDE_DIRS})
  set(WITH_HDF5 1)
else()
  set(HDF5_LIBRARIES)
  set(WITH_HDF5 0)
endif()

# find_package(Boost COMPONENTS regex)
if(NOT Boost_FOUND)
  find_package(Boost REQUIRED)
else()
  set(HAVE_BOOST_REGEX 1)
endif()

include_directories(${Boost_INCLUDE_DIRS})

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

include_directories(${CMAKE_BINARY_DIR})

# if (NOT NO_SHARED)
  add_subdirectory(greens_functions)  #XXX: submodule
# endif()

add_subdirectory(ecell4)
add_subdirectory(python)
