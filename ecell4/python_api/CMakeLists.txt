pybind11_add_module(ecell4 main.cpp)
target_link_libraries(ecell4 PRIVATE
    ecell4-core
    ecell4-bd
    ecell4-egfrd
    ecell4-gillespie
    ecell4-meso
    ecell4-ode
    ecell4-spatiocyte)

