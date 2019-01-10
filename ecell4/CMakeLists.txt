add_subdirectory(core)
add_subdirectory(egfrd)
add_subdirectory(gillespie)
add_subdirectory(ode)
add_subdirectory(bd)
add_subdirectory(meso)
add_subdirectory(spatiocyte)

pybind11_add_module(ecell4 python_api.cpp)
target_link_libraries(ecell4 PRIVATE
    ecell4-core
    ecell4-bd
    ecell4-egfrd
    ecell4-gillespie
    ecell4-meso
    ecell4-ode
    ecell4-spatiocyte)
