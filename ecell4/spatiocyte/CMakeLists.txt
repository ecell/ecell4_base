file(GLOB CPP_FILES *.cpp)

add_library(ecell4-spatiocyte STATIC ${CPP_FILES})
target_link_libraries(ecell4-spatiocyte INTERFACE ecell4-core)

add_subdirectory(tests)
add_subdirectory(samples)
