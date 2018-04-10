add_executable(hardbody hardbody.cpp)
target_link_libraries(hardbody ecell4-bd)

add_executable(gensphere generate_sphere_polygon.cpp)
target_link_libraries(gensphere ecell4-core ${GSL_LIBRARIES} ${GSL_CBLAS_LIBRARIES})

add_executable(surface surface.cpp)
target_link_libraries(surface ecell4-bd)
