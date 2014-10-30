set(CPP_FILES
    LatticeSimulator.cpp LatticeWorld.cpp)

set(HPP_FILES
    LatticeSimulator.hpp LatticeWorld.hpp LatticeFactory.hpp)

add_library(ecell4-lattice SHARED ${CPP_FILES} ${HPP_FILES})
target_link_libraries(ecell4-lattice ecell4-core)

add_subdirectory(tests)
add_subdirectory(samples)

install(TARGETS ecell4-lattice DESTINATION lib)
install(FILES ${HPP_FILES} DESTINATION "include/ecell4/lattice")
