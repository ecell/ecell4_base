if (NO_SHARED)
    return()
endif()

set(CPP_FILES
    GillespieSimulator.cpp GillespieWorld.cpp)

set(HPP_FILES
    GillespieSimulator.hpp GillespieWorld.hpp GillespieFactory.hpp)

add_library(ecell4-gillespie SHARED ${CPP_FILES} ${HPP_FILES})
target_link_libraries(ecell4-gillespie ecell4-core)

set(ECELL4_SHARED_DIRS ${CMAKE_CURRENT_BINARY_DIR}:${ECELL4_SHARED_DIRS} PARENT_SCOPE)

add_subdirectory(tests)
add_subdirectory(samples)

install(TARGETS ecell4-gillespie DESTINATION lib)
install(FILES ${HPP_FILES} DESTINATION "include/ecell4/gillespie")
