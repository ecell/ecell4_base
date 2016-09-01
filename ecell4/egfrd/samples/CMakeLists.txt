add_executable(mymapk mymapk.cpp)
target_link_libraries(mymapk ecell4-egfrd)

add_executable(polygon polygon.cpp)
target_link_libraries(polygon ecell4-egfrd)

add_executable(stl2vmd make_polygon_vmd.cpp)
target_link_libraries(stl2vmd ecell4-egfrd)
