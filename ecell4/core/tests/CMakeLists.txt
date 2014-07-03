set(TEST_NAMES
    Position3_test CompartmentSpace_test Species_test
    ReactionRule_test NetworkModel_test get_mapper_mf_test
    LatticeSpace_test EventScheduler_test)

foreach(TEST_NAME ${TEST_NAMES})
    add_executable(${TEST_NAME} ${TEST_NAME}.cpp)
    target_link_libraries(${TEST_NAME} ecell4-core boost_unit_test_framework.a)
    add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
endforeach(TEST_NAME)
