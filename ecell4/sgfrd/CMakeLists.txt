if (NO_SHARED)
    return()
endif()

file(GLOB HPP_FILES *.hpp)
set(CPP_FILES SGFRDWorld.cpp SGFRDSimulator.cpp)

# add_definitions(-DECELL4_SGFRD_NO_TRACE)

add_executable(sgfrd-sample-reaction sample/reaction.cpp ${CPP_FILES})
target_link_libraries(sgfrd-sample-reaction ecell4-egfrd ecell4-core ecell4-bd)
set_target_properties(sgfrd-sample-reaction PROPERTIES COMPILE_FLAGS "-std=c++11")

# add_library(ecell4-sgfrd SHARED ${CPP_FILES} ${HPP_FILES})
# target_link_libraries(ecell4-sgfrd ecell4-core ecell4-egfrd)
#
# set(ECELL4_SHARED_DIRS ${CMAKE_CURRENT_BINARY_DIR}:${ECELL4_SHARED_DIRS} PARENT_SCOPE)
#
# add_subdirectory(tests)
# add_subdirectory(samples)
#
# install(TARGETS ecell4-sgfrd DESTINATION lib)
# install(FILES ${HPP_FILES} DESTINATION "include/ecell4/sgfrd")
