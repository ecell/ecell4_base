add_executable(equilibrium equilibrium.cpp)
target_link_libraries(equilibrium ecell4-ode)

add_executable(dissociation dissociation.cpp)
target_link_libraries(dissociation ecell4-ode)

add_executable(equilibrium2 equilibrium2.cpp)
target_link_libraries(equilibrium2 ecell4-ode-ratelow)

add_executable(dissociation2 dissociation2.cpp)
target_link_libraries(dissociation2 ecell4-ode-ratelow)
