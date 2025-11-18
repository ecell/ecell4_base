cmake_minimum_required(VERSION 3.14)
project(ecell4 VERSION 4.0 LANGUAGES C CXX)

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
include_directories(${GSL_INCLUDE_DIRS})

# Robust pybind11 handling: use installed pybind11 if available, otherwise fetch it.
include(FetchContent)
find_package(pybind11 CONFIG QUIET)
if(NOT pybind11_FOUND)
  message(STATUS "pybind11 not found by find_package; fetching pybind11 with FetchContent")
  FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11.git
    # Pin a known stable tag; adjust as needed.
    GIT_TAG v2.11.0
  )
  FetchContent_MakeAvailable(pybind11)
endif()

# Only add subdirectories if they exist, to avoid configure-time failure for packaging contexts.
if(EXISTS "${CMAKE_SOURCE_DIR}/greens_functions/CMakeLists.txt")
  add_subdirectory(greens_functions)
else()
  message(STATUS "Skipping greens_functions (no CMakeLists.txt found)")
endif()

add_subdirectory(ecell4)
