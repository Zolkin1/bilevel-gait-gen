find_package(Catch2 3 REQUIRED)

add_executable(sim_test simulator_test.cpp)

message(STATUS "Making test directory.")

target_link_libraries(sim_test PRIVATE Catch2::Catch2WithMain)
target_link_libraries(sim_test PUBLIC simulator)
target_link_libraries(sim_test PUBLIC controller)
target_link_libraries(sim_test PUBLIC utils)

add_executable(mpc_test mpc_test.cpp)
target_link_libraries(mpc_test PRIVATE Catch2::Catch2WithMain)
target_link_libraries(mpc_test PUBLIC mpc)
target_link_libraries(mpc_test PUBLIC utils)

add_executable(spline_test splines_tests.cpp)
target_link_libraries(spline_test PRIVATE Catch2::Catch2WithMain)
target_link_libraries(spline_test PUBLIC mpc)
target_link_libraries(spline_test PUBLIC utils)

include(Catch)
catch_discover_tests(sim_test)
catch_discover_tests(mpc_test)

add_custom_command(TARGET sim_test POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${CMAKE_CURRENT_SOURCE_DIR}/a1_configuration_test.yaml
        ${CMAKE_CURRENT_SOURCE_DIR}/test_config.yaml
        $<TARGET_FILE_DIR:sim_test>
)

add_custom_command(TARGET mpc_test POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${CMAKE_CURRENT_SOURCE_DIR}/a1_configuration_test.yaml
        $<TARGET_FILE_DIR:sim_test>
)

add_executable(gait_opt_playground gait_opt_playground.cpp)
target_link_libraries(gait_opt_playground utils)
target_link_libraries(gait_opt_playground mpc)
target_link_libraries(gait_opt_playground simulator)

add_executable(gait_opt_ls gait_opt_line_search.cpp)
target_link_libraries(gait_opt_ls utils)
target_link_libraries(gait_opt_ls mpc)
target_link_libraries(gait_opt_ls simulator)

add_executable(sim_mpc simulation_mpc.cpp)
target_link_libraries(sim_mpc utils)
target_link_libraries(sim_mpc mpc)
target_link_libraries(sim_mpc simulator)

#catch_discover_tests(gait_opt_playground.cpp)