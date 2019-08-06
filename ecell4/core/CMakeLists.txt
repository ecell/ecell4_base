set(ECELL4_VERSION ${PROJECT_VERSION})

configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/config.h.in
    ${CMAKE_CURRENT_BINARY_DIR}/config.h)

# file(READ ${CMAKE_CURRENT_BINARY_DIR}/config.h ECELL4_GENERATED_CONFIG_H)
# message(STATUS ${ECELL4_GENERATED_CONFIG_H})

file(GLOB CPP_FILES *.cpp)

add_library(ecell4-core STATIC ${CPP_FILES})

target_link_libraries(ecell4-core PRIVATE
    ${HDF5_LIBRARIES} ${Boost_LIBRARIES} ${GSL_LIBRARIES} ${GSL_CBLAS_LIBRARIES})

if(WITH_VTK AND NOT VTK_LIBRARIES)
    target_link_libraries(ecell4-core PRIVATE vtkHybrid vtkWidgets)
else()
    target_link_libraries(ecell4-core PRIVATE ${VTK_LIBRARIES})
endif()

add_subdirectory(tests)
