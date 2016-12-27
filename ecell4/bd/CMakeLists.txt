if (NO_SHARED)
    return()
endif()

set(CPP_FILES
    BDSimulator.cpp BDPropagator.cpp functions3d.cpp
    BDContainer2D.cpp BDPolygon.cpp BDPropagator2D.cpp rotate_vector.cpp
    )

set(HPP_FILES
    BDSimulator.hpp BDWorld.hpp BDPropagator.hpp functions3d.hpp BDFactory.hpp
    BDContainer2D.hpp BDPolygon.hpp BDPropagator2D.hpp rotate_vector.hpp
    )

add_library(ecell4-bd SHARED ${CPP_FILES} ${HPP_FILES})
target_link_libraries(ecell4-bd ecell4-core)

set(ECELL4_SHARED_DIRS ${CMAKE_CURRENT_BINARY_DIR}:${ECELL4_SHARED_DIRS} PARENT_SCOPE)

add_subdirectory(tests)
add_subdirectory(samples)

install(TARGETS ecell4-bd DESTINATION lib)
install(FILES ${HPP_FILES} DESTINATION "include/ecell4/bd")
