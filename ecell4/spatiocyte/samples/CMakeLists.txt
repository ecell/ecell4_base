add_executable(simple_lattice simple_lattice.cpp)
target_link_libraries(simple_lattice ecell4-spatiocyte)

add_executable(diffusion diffusion.cpp)
target_link_libraries(diffusion ecell4-spatiocyte)
