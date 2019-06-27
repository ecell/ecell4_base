file(GLOB CPP_FILES *.cpp)

add_library(ecell4-sgfrd STATIC ${CPP_FILES})
target_compile_definitions(ecell4-sgfrd PUBLIC -DECELL4_SGFRD_NO_TRACE)
target_link_libraries(ecell4-sgfrd INTERFACE ecell4-core greens_functions)

add_subdirectory(samples)
