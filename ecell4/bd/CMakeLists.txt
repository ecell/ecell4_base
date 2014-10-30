set(CPP_FILES
    BDSimulator.cpp BDPropagator.cpp functions3d.cpp)

set(HPP_FILES
    BDSimulator.hpp BDWorld.hpp BDPropagator.hpp functions3d.hpp BDFactory.hpp)

add_library(ecell4-bd SHARED ${CPP_FILES} ${HPP_FILES})
target_link_libraries(ecell4-bd ecell4-core)

add_subdirectory(tests)
add_subdirectory(samples)

install(TARGETS ecell4-bd DESTINATION lib)
install(FILES ${HPP_FILES} DESTINATION "include/ecell4/bd")
