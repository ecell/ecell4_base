add_subdirectory(greens_functions)

add_subdirectory(core)
add_subdirectory(egfrd)
add_subdirectory(gillespie)
add_subdirectory(ode)
add_subdirectory(bd)
add_subdirectory(meso)
add_subdirectory(spatiocyte)

if (NOT NO_SHARED)
  set(ECELL4_SHARED_DIRS ${ECELL4_SHARED_DIRS} PARENT_SCOPE)
endif()
