# SET(CPACK_GENERATOR "DEB")
# SET(CPACK_DEBIAN_PACKAGE_DEPENDS "libgsl0-dev (>= 1.16+dfsg-1ubuntu1), libhdf5-serial-dev (>= 1.8.11-5ubuntu7), libboost-dev (>= 1.54.0.1ubuntu1)")
# SET(CPACK_PACKAGE_CONTACT "Kozo Nishida <knishida@riken.jp>")
# SET(CPACK_DEBIAN_PACKAGE_MAINTAINER "Kozo Nishida <knishida@riken.jp>")
# INCLUDE(CPack)

configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/config.h.in
    ${CMAKE_CURRENT_BINARY_DIR}/config.h)

if (NO_SHARED)
    return()
endif()

file(GLOB CPP_FILES *.cpp)
file(GLOB HPP_FILES *.hpp)

add_library(ecell4-core SHARED
    ${CPP_FILES} ${HPP_FILES} "${CMAKE_CURRENT_BINARY_DIR}/config.h")

set(ECELL4_SHARED_DIRS ${CMAKE_CURRENT_BINARY_DIR}:${ECELL4_SHARED_DIRS} PARENT_SCOPE)

if(WITH_VTK AND NOT VTK_LIBRARIES)
    target_link_libraries(ecell4-core
        ${HDF5_LIBRARIES} ${Boost_LIBRARIES} ${GSL_LIBRARIES} ${GSL_CBLAS_LIBRARIES}
        vtkHybrid vtkWidgets)
else()
    target_link_libraries(ecell4-core
        ${HDF5_LIBRARIES} ${Boost_LIBRARIES} ${GSL_LIBRARIES} ${GSL_CBLAS_LIBRARIES}
        ${VTK_LIBRARIES})
endif()

add_subdirectory(tests)

install(TARGETS ecell4-core DESTINATION lib)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/config.h" ${HPP_FILES}
    DESTINATION "include/ecell4/core")
