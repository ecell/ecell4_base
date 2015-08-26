set(CPP_FILES
	#ODESimulator.cpp 
	ODEWorld.cpp ODENetworkModel.cpp ODEReactionRule.cpp ODESimulator2.cpp ODERatelaw.cpp)

set(HPP_FILES
	#ODESimulator.hpp 
	ODEWorld.hpp ODEFactory.hpp ODEReactionRule.hpp ODERatelaw.hpp ODENetworkModel.hpp ODESimulator2.hpp ODEFactory2.hpp)

add_library(ecell4-ode SHARED ${CPP_FILES} ${HPP_FILES})
target_link_libraries(ecell4-ode ecell4-core)

#add_subdirectory(tests)
#add_subdirectory(samples)

install(TARGETS ecell4-ode DESTINATION lib)
install(FILES ${HPP_FILES} DESTINATION "include/ecell4/ode")
