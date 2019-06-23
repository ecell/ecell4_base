pybind11_add_module(ecell4_base
    main.cpp
    core.cpp
    bd.cpp
    egfrd.cpp
    gillespie.cpp
    meso.cpp
    ode.cpp
    sgfrd.cpp
    spatiocyte.cpp)
target_link_libraries(ecell4_base PRIVATE
    ecell4-core
    ecell4-bd
    ecell4-egfrd
    ecell4-gillespie
    ecell4-meso
    ecell4-ode
    ecell4-sgfrd
    ecell4-spatiocyte)

set_target_properties(ecell4_base PROPERTIES BUILD_WITH_INSTALL_RPATH TRUE)

if (UNIX AND NOT APPLE)
    set_target_properties(ecell4_base PROPERTIES INSTALL_RPATH "$ORIGIN")
elseif (APPLE)
    set_target_properties(ecell4_base PROPERTIES INSTALL_RPATH "@loader_path")
endif ()
