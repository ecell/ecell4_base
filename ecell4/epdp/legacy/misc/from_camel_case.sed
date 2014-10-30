# ls *.py | xargs sed -i -r -f misc/from_camel_case.sed
# Don't run this in gui/, too many vtk function calls.

s/(\<[^A-Z][a-z_]+)([A-Z][a-z]+)/\1_\l\2/g
s/(\<[^A-Z][a-z_]+)([A-Z][a-z]+)/\1_\l\2/g
s/(\<[^A-Z][a-z_]+)([A-Z][a-z]+)/\1_\l\2/g
s/(\<[^A-Z][a-z_]+)([A-Z][a-z]+)/\1_\l\2/g
s/(\<[^A-Z][a-z_]+)([A-Z][a-z]+)/\1_\l\2/g
s/(\<[^A-Z][a-z_]+)([A-Z][a-z]+)/\1_\l\2/g
s/(\<[^A-Z][a-z_]+)([A-Z][a-z]+)/\1_\l\2/g


# Extra.

s/CoM/com/g
s/coM/com/g
s/vectorS/vector_s/g
s/vectorX/vector_x/g
s/vectorY/vector_y/g
s/ZAxis/_z_axis/g
s/unitX/unit_x/g
s/unitY/unit_y/g
s/unitZ/unit_z/g
s/domainIDGenerator/domain_id_generator/g
s/shellIDGenerator/shell_id_generator/g
s/drawBDdisplacement/draw_bd_displacement/g
s/calculateBDDt/calculate_bd_dt/g
s/jTable/j_table/g
s/yTable/y_table/g
s/zTable/z_table/g
s/new_inter_particleS/new_inter_particle_s/g
s/newIV/new_iv/g
s/eventID/event_id/g

#s/particleA/particle_a/g
#s/particleB/particle_b/g
#s/shell_id0/shell_id_0/g
#s/shell_id1/shell_id_1/g
#s/test_alphan/test_alpha_n/g

#s/tR/t_R/g
#s/aR/a_R/g

#s/kD/k_D/g
#s/getD/get_D/g
#s/drawR/draw_r/g


# Revert

# Python Logger.
s/get_logger/getLogger/g
s/set_formatter/setFormatter/g
s/add_handler/addHandler/g
s/set_level/setLevel/g
s/max_bytes/maxBytes/g

# Python UnitTest.
s/set_up/setUp/g
s/setUp_base/setUpBase/g
s/tear_down/tearDown/g
s/fail_if/failIf/g
s/test__/test_/g
s/assert_equal/assertEqual/g
s/assert_not_equal/assertNotEqual/g
s/assert_almost_equal/assertAlmostEqual/g
s/assert_not_almost_equal/assertNotAlmostEqual/g
s/assert_true/assertTrue/g
s/assert_raises/assertRaises/g

# Samples variables.
s/k_f__rp/k_f_Rp/g
s/k_b__rp/k_b_Rp/g
s/fraction__sp/fraction_Sp/g

# Ecell.
s/get_data/getData/g
s/load_model/loadModel/g
s/create_logger_stub/createLoggerStub/g
s/create_entity_stub/createEntityStub/g
s/get_stdout/getStdout/g
s/register_ecell_session/registerEcellSession/g

# EventScheduler.
s/add_event/addEvent/g
s/get_top_event/getTopEvent/g
s/get_time/getTime/g
s/get_arg/getArg/g
s/get_top_time/getTopTime/g
s/get_topID/getTopID/g
s/get_size/getSize/g
s/get_event_by_index/getEventByIndex/g
s/update_event_time/updateEventTime/g
s/remove_event/removeEvent/g

# World.
s/calculate_pair__com/calculate_pair_CoM/g

# Green's functions.
s/draw_time/drawTime/g
s/draw_event_type/drawEventType/g
s/draw_theta/drawTheta/g
s/get_sigma/getSigma/g

# Other.
s/array__simple/array_simple/g
s/array__cyclic/array_cyclic/g
s/bessel__ynu/bessel_Ynu/g

# Extra after revert
s/test_drawTime/test_draw_time/g
s/test_drawEventType/test_draw_event_type/g
s/test_drawTheta/test_draw_theta/g
#s/test_drawR/test_draw_r/g
