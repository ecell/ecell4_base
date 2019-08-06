cmake_minimum_required(VERSION 3.1)

project(ecell4)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -pg -O0 -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2")

if(WIN32)
    add_definitions(-DNOMINMAX)
    add_definitions(-DGSL_DLL)
endif()

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

message(STATUS "Looking for HDF5")
find_package(HDF5 1.10 COMPONENTS C CXX HL REQUIRED)
if (HDF5_FOUND)
  message(STATUS "HDF5 was found ... ${HDF5_VERSION}")
  include_directories(${HDF5_INCLUDE_DIRS})
  set(WITH_HDF5 1)
else()
  message(STATUS "Could NOT find HDF5")
  set(HDF5_LIBRARIES)
  set(WITH_HDF5 0)
endif()

find_package(Boost)
if(NOT DEFINED Boost_VERSION)
    message(FATAL_ERROR
        "ecell4 requires Boost C++ Library. "
        "If you have already installed it, try `-DBOOST_ROOT=/path/to/boost`")
else()
    include_directories(${Boost_INCLUDE_DIRS})
    link_directories(${Boost_LIBRARY_DIRS})
endif()

find_package(GSL REQUIRED)
include_directories({${GSL_INCLUDE_DIRS})

include(CheckIncludeFileCXX)
include(CheckCXXSourceCompiles)

check_include_file_cxx(unordered_map HAVE_UNORDERED_MAP)
check_cxx_source_compiles("#include <tr1/unordered_map>
int main() { std::tr1::unordered_map<int, int> a, b(a); return 0; }"
    HAVE_TR1_UNORDERED_MAP)
check_include_file_cxx("boost/unordered_map.hpp" HAVE_BOOST_UNORDERED_MAP_HPP)

check_cxx_source_compiles("#include <tr1/functional>
int main() { std::tr1::hash<int>(); return 0; }"
    HAVE_TR1_FUNCTIONAL)
check_cxx_source_compiles("#include <functional>
int main() { std::hash<int>(); return 0; }"
    HAVE_STD_HASH)
check_include_file_cxx("boost/functional/hash.hpp" HAVE_BOOST_FUNCTIONAL_HASH_HPP)

check_include_file_cxx(chrono HAVE_CHRONO)

add_subdirectory(greens_functions)
add_subdirectory(pybind11)
add_subdirectory(ecell4)
