add_executable(hardbody hardbody.cpp)
target_link_libraries(hardbody ecell4-bd)

add_executable(surface surface.cpp)
target_link_libraries(surface ecell4-bd)
