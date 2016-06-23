# cmake_minimum_required(VERSION 2.8)
# INCLUDE(CheckCXXSourceCompiles)

add_executable(make_sjy_table make_sjy_table.cpp)
add_executable(make_cjy_table make_cjy_table.cpp)
target_link_libraries(make_sjy_table gslcblas)
target_link_libraries(make_cjy_table gsl gslcblas)
# add_custom_command(
#     OUTPUT SphericalBesselTable.hpp
#     COMMAND make_sjy_table ${PROJECT_BINARY_DIR}/SphericalBesselTable.hpp
#     DEPENDS make_sjy_table)
# add_custom_command(
#     OUTPUT CylindricalBesselTable.hpp
#     COMMAND make_cjy_table ${PROJECT_BINARY_DIR}/CylindricalBesselTable.hpp
#     DEPENDS make_cjy_table)
# add_custom_target(BesselTables ALL DEPENDS ${PROJECT_BINARY_DIR}/SphericalBesselTable.hpp ${PROJECT_BINARY_DIR}/CylindricalBesselTable.hpp)
add_custom_command(
    OUTPUT ${PROJECT_BINARY_DIR}/ecell4/egfrd/SphericalBesselTable.hpp
    COMMAND make_sjy_table ${PROJECT_BINARY_DIR}/ecell4/egfrd/SphericalBesselTable.hpp
    DEPENDS make_sjy_table)
add_custom_command(
    OUTPUT ${PROJECT_BINARY_DIR}/ecell4/egfrd/CylindricalBesselTable.hpp
    COMMAND make_cjy_table ${PROJECT_BINARY_DIR}/ecell4/egfrd/CylindricalBesselTable.hpp
    DEPENDS make_cjy_table)
add_custom_target(BesselTables ALL DEPENDS ${PROJECT_BINARY_DIR}/ecell4/egfrd/SphericalBesselTable.hpp ${PROJECT_BINARY_DIR}/ecell4/egfrd/CylindricalBesselTable.hpp)
