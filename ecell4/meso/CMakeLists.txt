set(CPP_FILES
    MesoscopicSimulator.cpp MesoscopicWorld.cpp)

set(HPP_FILES
    MesoscopicSimulator.hpp MesoscopicWorld.hpp)

add_library(ecell4-meso SHARED ${CPP_FILES} ${HPP_FILES})
target_link_libraries(ecell4-meso ecell4-core)

add_subdirectory(tests)
# add_subdirectory(samples)

install(TARGETS ecell4-meso DESTINATION lib)
install(FILES ${HPP_FILES} DESTINATION "include/ecell4/meso")
