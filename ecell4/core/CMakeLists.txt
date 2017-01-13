# SET(CPACK_GENERATOR "DEB")
# SET(CPACK_DEBIAN_PACKAGE_DEPENDS "libgsl0-dev (>= 1.16+dfsg-1ubuntu1), libhdf5-serial-dev (>= 1.8.11-5ubuntu7), libboost-dev (>= 1.54.0.1ubuntu1)")
# SET(CPACK_PACKAGE_CONTACT "Kozo Nishida <knishida@riken.jp>")
# SET(CPACK_DEBIAN_PACKAGE_MAINTAINER "Kozo Nishida <knishida@riken.jp>")
# INCLUDE(CPack)

set(ECELL4_VERSION ${PROJECT_VERSION})

configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/config.h.in
    ${CMAKE_CURRENT_BINARY_DIR}/config.h)

if (NO_SHARED)
    return()
endif()

# set(CPP_FILES
#     CompartmentSpace.cpp Context.cpp Integer3.cpp Journal.cpp ReactionRule.cpp
#     LatticeSpace.cpp Model.cpp NetworkModel.cpp NetfreeModel.cpp ParticleSpace.cpp
#     Real3.cpp RandomNumberGenerator.cpp Species.cpp UnitSpecies.cpp
#     Sphere.cpp Cylinder.cpp SubvolumeSpace.cpp ParticleSpaceCellListImpl.cpp
#     LatticeSpaceCellListImpl.cpp collision.cpp Rod.cpp AABB.cpp AABBSurface.cpp
#     PlanarSurface.cpp Mesh.cpp observers.cpp extras.cpp)

# set(HPP_FILES
#     CompartmentSpace.hpp DynamicPriorityQueue.hpp LatticeSpace.hpp
#     MolecularType.hpp Particle.hpp Real3.hpp ReactionRule.hpp Space.hpp
#     VacantType.hpp extras.hpp linear_algebra.hpp CompartmentSpaceHDF5Writer.hpp
#     EventScheduler.hpp Identifier.hpp LatticeSpaceHDF5Writer.hpp
#     SubvolumeSpace.hpp SubvolumeSpaceHDF5Writer.hpp
#     MolecularTypeBase.hpp ParticleSpace.hpp RandomNumberGenerator.hpp
#     SerialIDGenerator.hpp Species.hpp Voxel.hpp functions.hpp swap.hpp ModelWrapper.hpp
#     Context.hpp Integer3.hpp Journal.hpp Model.hpp NetworkModel.hpp NetfreeModel.hpp
#     ParticleSpaceHDF5Writer.hpp Reaction.hpp Simulator.hpp SimulatorBase.hpp
#     UnitSpecies.hpp observers.hpp exceptions.hpp get_mapper_mf.hpp types.hpp
#     Shape.hpp Sphere.hpp Cylinder.hpp AABB.hpp AABBSurface.hpp StructureType.hpp
#     ParticleSpaceCellListImpl.hpp LatticeSpaceCellListImpl.hpp
#     PlanarSurface.hpp Rod.hpp collision.hpp
#     comparators.hpp SimulatorFactory.hpp Mesh.hpp shape_operators.hpp)

file(GLOB CPP_FILES *.cpp)
file(GLOB HPP_FILES *.hpp)

add_library(ecell4-core SHARED
    ${CPP_FILES} ${HPP_FILES} "${CMAKE_CURRENT_BINARY_DIR}/config.h")

set(ECELL4_SHARED_DIRS ${CMAKE_CURRENT_BINARY_DIR}:${ECELL4_SHARED_DIRS} PARENT_SCOPE)

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
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/config.h" ${HPP_FILES}
    DESTINATION "include/ecell4/core")
