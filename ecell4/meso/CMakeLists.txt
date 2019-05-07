file(GLOB CPP_FILES *.cpp)

add_library(ecell4-meso STATIC ${CPP_FILES})
target_link_libraries(ecell4-meso INTERFACE ecell4-core)

add_subdirectory(tests)
add_subdirectory(samples)
