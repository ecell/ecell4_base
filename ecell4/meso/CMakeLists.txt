if (NO_SHARED)
    return()
endif()

set(CPP_FILES
    MesoscopicSimulator.cpp MesoscopicWorld.cpp)

set(HPP_FILES
    MesoscopicSimulator.hpp MesoscopicWorld.hpp MesoscopicFactory.hpp)

add_library(ecell4-meso SHARED ${CPP_FILES} ${HPP_FILES})
target_link_libraries(ecell4-meso ecell4-core)

set(ECELL4_SHARED_DIRS ${CMAKE_CURRENT_BINARY_DIR}:${ECELL4_SHARED_DIRS} PARENT_SCOPE)

add_subdirectory(tests)
add_subdirectory(samples)

install(TARGETS ecell4-meso DESTINATION lib)
install(FILES ${HPP_FILES} DESTINATION "include/ecell4/meso")
