set(CPP_FILES
    SpatiocyteSimulator.cpp SpatiocyteWorld.cpp)

set(HPP_FILES
    SpatiocyteSimulator.hpp SpatiocyteWorld.hpp SpatiocyteFactory.hpp)

set(BINARY_DIRECTORIES ${CMAKE_CURRENT_BINARY_DIR}:${BINARY_DIRECTORIES} PARENT_SCOPE)

add_library(ecell4-spatiocyte SHARED ${CPP_FILES} ${HPP_FILES})
target_link_libraries(ecell4-spatiocyte ecell4-core)

add_subdirectory(tests)
add_subdirectory(samples)

install(TARGETS ecell4-spatiocyte DESTINATION lib)
install(FILES ${HPP_FILES} DESTINATION "include/ecell4/spatiocyte")
