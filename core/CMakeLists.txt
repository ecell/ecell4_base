cmake_minimum_required(VERSION 2.6)
project(ecell4-core)

set(CPP_FILES
    CompartmentSpace.cpp Context.cpp Global.cpp Journal.cpp
    LatticeSpace.cpp Model.cpp NetworkModel.cpp ParticleSpace.cpp
    Position3.cpp RandomNumberGenerator.cpp Species.cpp UnitSpecies.cpp)

include(CheckIncludeFiles)
check_include_files(unordered_map.hpp HAVE_UNORDERED_MAP_HPP)
check_include_files("boost/unordered_map.hpp" HAVE_BOOST_UNORDERED_MAP_HPP)
check_include_files("boost/functional/hash.hpp" HAVE_BOOST_FUNCTIONAL_HASH_HPP)

include(CheckCXXSourceCompiles)
check_cxx_source_compiles("#include <tr1/unordered_map>
int main() { std::tr1::unordered_map<int, int> a, b(a); return 0; }"
    HAVE_TR1_UNORDERED_MAP)
check_cxx_source_compiles("#include <tr1/functional>
int main() { std::tr1::hash<int>(); return 0; }"
    HAVE_TR1_FUNCTIONAL)

find_package(HDF5 COMPONENTS C CXX HL REQUIRED)
include_directories(${HDF5_INCLUDE_DIRS})
find_package(Boost COMPONENTS regex REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

find_library(GSL_LIBRARIES gsl)
find_library(GSL_CBLAS_LIBRARIES gslcblas)

configure_file(
    ${PROJECT_SOURCE_DIR}/config.h.in
    ${PROJECT_BINARY_DIR}/config.h)
include_directories(${PROJECT_BINARY_DIR})

add_library(ecell4-core SHARED ${CPP_FILES})
# target_link_libraries(ecell4-core ${HDF5_LIBRARIES} ${Boost_LIBRARIES} gsl gslcblas)
target_link_libraries(ecell4-core ${HDF5_LIBRARIES} ${Boost_LIBRARIES} ${GSL_LIBRARIES} ${GSL_CBLAS_LIBRARIES})

enable_testing()

add_subdirectory(tests)
