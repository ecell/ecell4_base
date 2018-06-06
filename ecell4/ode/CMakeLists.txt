if (NO_SHARED)
    return()
endif()

set(CPP_FILES
	ODEWorld.cpp ODENetworkModel.cpp ODEReactionRule.cpp ODESimulator.cpp ODERatelaw.cpp ODESimulator_New.cpp ODEWorld_New.cpp)

set(HPP_FILES
	ODEWorld.hpp ODEFactory.hpp ODEReactionRule.hpp ODERatelaw.hpp ODENetworkModel.hpp ODESimulator.hpp ODEFactory.hpp ODESimulator_New.hpp ODEWorld_New.hpp)

add_library(ecell4-ode SHARED ${CPP_FILES} ${HPP_FILES})
target_link_libraries(ecell4-ode ecell4-core)

set(ECELL4_SHARED_DIRS ${CMAKE_CURRENT_BINARY_DIR}:${ECELL4_SHARED_DIRS} PARENT_SCOPE)

add_subdirectory(tests)
add_subdirectory(samples)

install(TARGETS ecell4-ode DESTINATION lib)
install(FILES ${HPP_FILES} DESTINATION "include/ecell4/ode")
