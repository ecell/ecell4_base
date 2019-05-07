file(GLOB CPP_FILES *.cpp)
add_definitions(-DECELL4_SGFRD_NO_TRACE)

add_library(ecell4-sgfrd STATIC ${CPP_FILES})
target_link_libraries(ecell4-sgfrd INTERFACE ecell4-core)

add_subdirectory(samples)
