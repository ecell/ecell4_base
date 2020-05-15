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
  if(WIN32)
    add_definitions(-DH5_BUILT_AS_DYNAMIC_LIB)
  endif()
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

find_package(pybind11)
if(NOT pybind11_FOUND)
    add_subdirectory(pybind11)
endif()

add_subdirectory(greens_functions)
add_subdirectory(ecell4)
