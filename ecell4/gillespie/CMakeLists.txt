set(CPP_FILES
    GillespieSimulator.cpp GillespieWorld.cpp)

set(HPP_FILES
    GillespieSimulator.hpp GillespieWorld.hpp GillespieFactory.hpp)

add_library(ecell4-gillespie SHARED ${CPP_FILES} ${HPP_FILES})
target_link_libraries(ecell4-gillespie ecell4-core)

add_subdirectory(tests)
add_subdirectory(samples)

install(TARGETS ecell4-gillespie DESTINATION lib)
install(FILES ${HPP_FILES} DESTINATION "include/ecell4/gillespie")
