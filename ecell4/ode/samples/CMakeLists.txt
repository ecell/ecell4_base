add_executable(equilibrium equilibrium.cpp)
target_link_libraries(equilibrium ecell4-ode)

add_executable(dissociation dissociation.cpp)
target_link_libraries(dissociation ecell4-ode)

