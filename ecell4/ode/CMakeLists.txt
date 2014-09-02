set(CPP_FILES
    ODESimulator.cpp ODEWorld.cpp)

set(HPP_FILES
    ODESimulator.hpp ODEWorld.hpp)

set(RATELOW_CPP_FILES
	ODESimulator_ratelow.cpp ODEWorld.cpp)
set(RATELOW_HPP_FILES
	ODESimulator_ratelow.hpp ODEWorld.hpp)

add_library(ecell4-ode SHARED ${CPP_FILES} ${HPP_FILES})
target_link_libraries(ecell4-ode ecell4-core)

add_library(ecell4-ode-ratelow SHARED ${RATELOW_CPP_FILES} ${RATELOW_HPP_FILES})
target_link_libraries(ecell4-ode-ratelow ecell4-core)

add_subdirectory(tests)
add_subdirectory(samples)

install(TARGETS ecell4-ode DESTINATION lib)
install(FILES ${HPP_FILES} DESTINATION "include/ecell4/ode")
