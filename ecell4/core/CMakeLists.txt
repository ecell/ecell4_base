# SET(CPACK_GENERATOR "DEB")
# SET(CPACK_DEBIAN_PACKAGE_DEPENDS "libgsl0-dev (>= 1.16+dfsg-1ubuntu1), libhdf5-serial-dev (>= 1.8.11-5ubuntu7), libboost-dev (>= 1.54.0.1ubuntu1)")
# SET(CPACK_PACKAGE_CONTACT "Kozo Nishida <knishida@riken.jp>")
# SET(CPACK_DEBIAN_PACKAGE_MAINTAINER "Kozo Nishida <knishida@riken.jp>")
# INCLUDE(CPack)

set(CPP_FILES
    CompartmentSpace.cpp Context.cpp Integer3.cpp Journal.cpp ReactionRule.cpp
    LatticeSpace.cpp Model.cpp NetworkModel.cpp NetfreeModel.cpp ParticleSpace.cpp
    Real3.cpp RandomNumberGenerator.cpp Species.cpp UnitSpecies.cpp
    Sphere.cpp Cylinder.cpp SubvolumeSpace.cpp ParticleSpaceCellListImpl.cpp
    LatticeSpaceCellListImpl.cpp collision.cpp Rod.cpp AABB.cpp PlanarSurface.cpp
    Mesh.cpp observers.cpp extras.cpp)

set(HPP_FILES
    CompartmentSpace.hpp DynamicPriorityQueue.hpp LatticeSpace.hpp
    MolecularType.hpp Particle.hpp Real3.hpp ReactionRule.hpp Space.hpp 
    VacantType.hpp extras.hpp linear_algebra.hpp CompartmentSpaceHDF5Writer.hpp
    EventScheduler.hpp Identifier.hpp LatticeSpaceHDF5Writer.hpp
    SubvolumeSpace.hpp SubvolumeSpaceHDF5Writer.hpp
    MolecularTypeBase.hpp ParticleSpace.hpp RandomNumberGenerator.hpp
    SerialIDGenerator.hpp Species.hpp Voxel.hpp functions.hpp swap.hpp ModelWrapper.hpp
    Context.hpp Integer3.hpp Journal.hpp Model.hpp NetworkModel.hpp NetfreeModel.hpp
    ParticleSpaceHDF5Writer.hpp Reaction.hpp Simulator.hpp SimulatorBase.hpp
    UnitSpecies.hpp observers.hpp exceptions.hpp get_mapper_mf.hpp types.hpp
    Shape.hpp Sphere.hpp Cylinder.hpp AABB.hpp StructureType.hpp InterfaceType.hpp
    ParticleSpaceCellListImpl.hpp LatticeSpaceCellListImpl.hpp
    PlanarSurface.hpp Rod.hpp collision.hpp
    comparators.hpp SimulatorFactory.hpp Mesh.hpp shape_operators.hpp)

if(WITH_VTK)
  find_package(VTK)
  if(VTK_FOUND)
    include(${VTK_USE_FILE})
    set(HAVE_VTK 1)
  endif()
endif()

if(WITH_HDF5)
find_package(HDF5 COMPONENTS C CXX HL REQUIRED)
include_directories(${HDF5_INCLUDE_DIRS})
set(WITH_HDF5 1)
endif()

# find_package(Boost COMPONENTS regex)
find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
# set(HAVE_BOOST_REGEX ${Boost_FOUND})

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

configure_file(
    ${PROJECT_SOURCE_DIR}/ecell4/core/config.h.in
    ${PROJECT_BINARY_DIR}/ecell4/core/config.h)

add_library(ecell4-core SHARED
    ${CPP_FILES} ${HPP_FILES} "${PROJECT_BINARY_DIR}/ecell4/core/config.h")

if(WITH_VTK AND NOT VTK_LIBRARIES)
    target_link_libraries(ecell4-core
        ${HDF5_LIBRARIES} ${Boost_LIBRARIES} ${GSL_LIBRARIES} ${GSL_CBLAS_LIBRARIES}
        vtkHybrid vtkWidgets)
else()
    target_link_libraries(ecell4-core
        ${HDF5_LIBRARIES} ${Boost_LIBRARIES} ${GSL_LIBRARIES} ${GSL_CBLAS_LIBRARIES}
        ${VTK_LIBRARIES})
endif()

add_subdirectory(tests)

install(TARGETS ecell4-core DESTINATION lib)
install(FILES "${PROJECT_BINARY_DIR}/ecell4/core/config.h" ${HPP_FILES}
    DESTINATION "include/ecell4/core")
