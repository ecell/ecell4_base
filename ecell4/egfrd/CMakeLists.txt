file(GLOB CPP_FILES *.cpp)

add_library(ecell4-egfrd STATIC ${CPP_FILES})
target_link_libraries(ecell4-egfrd INTERFACE ecell4-core)
target_link_libraries(ecell4-egfrd PUBLIC ${GSL_LIBRARIES} ${GSL_CBLAS_LIBRARIES} greens_functions)

# add_subdirectory(tests)
add_subdirectory(samples)
