pybind11_add_module(ecell4
    main.cpp
    core.cpp
    bd.cpp
    egfrd.cpp
    gillespie.cpp
    spatiocyte.cpp)
target_link_libraries(ecell4 PRIVATE
    ecell4-core
    ecell4-bd
    ecell4-egfrd
    ecell4-gillespie
    ecell4-meso
    ecell4-ode
    ecell4-spatiocyte)

set_target_properties(ecell4 PROPERTIES BUILD_WITH_INSTALL_RPATH TRUE)

if (UNIX AND NOT APPLE)
    set_target_properties(ecell4 PROPERTIES INSTALL_RPATH "$ORIGIN")
elseif (APPLE)
    set_target_properties(ecell4 PROPERTIES INSTALL_RPATH "@loader_path")
endif ()
