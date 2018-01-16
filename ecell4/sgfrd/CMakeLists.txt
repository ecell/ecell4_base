if (NO_SHARED)
    return()
endif()

file(GLOB HPP_FILES *.hpp)
set(CPP_FILES polygon_traits.cpp SGFRDWorld.cpp SGFRDSimulator.cpp)

# add_definitions(-DECELL4_SGFRD_NO_TRACE)

if( Boost_LIBRARIES )
    set(SGFRD_LIBRARIES_TO_BE_LINKED ${Boost_LIBRARIES})
else( Boost_LIBRARIES )
    # sometimes ${Boost_LIBRARIES} become empty so here make Boost.Chrono and
    # Boost.System header-only.
    add_definitions(-DBOOST_CHRONO_HEADER_ONLY)
    add_definitions(-DBOOST_SYSTEM_NO_DEPRECATED)
    add_definitions(-DBOOST_ERROR_CODE_HEADER_ONLY)
    add_definitions(-DBOOST_SYSTEM_NO_LIB)
endif( Boost_LIBRARIES )

add_executable(sgfrd-sample-diffusion sample/diffusion.cpp ${CPP_FILES})
target_link_libraries(sgfrd-sample-diffusion ecell4-egfrd ecell4-core ecell4-bd ${SGFRD_LIBRARIES_TO_BE_LINKED})

add_executable(sgfrd-sample-reaction sample/reaction.cpp ${CPP_FILES})
target_link_libraries(sgfrd-sample-reaction ecell4-egfrd ecell4-core ecell4-bd ${SGFRD_LIBRARIES_TO_BE_LINKED})

add_executable(sgfrd-sample-bdonly sample/bdonly.cpp ${CPP_FILES})
target_link_libraries(sgfrd-sample-bdonly ecell4-egfrd ecell4-core ecell4-bd ${SGFRD_LIBRARIES_TO_BE_LINKED})

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
