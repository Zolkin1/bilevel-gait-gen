find_package(Catch2 3 REQUIRED)

add_executable(controller_sim_test controller_test.cpp)

target_link_libraries(controller_sim_test PRIVATE Catch2::Catch2WithMain)

target_link_libraries(controller_sim_test PUBLIC simulator)
target_link_libraries(controller_sim_test PUBLIC mpc)

include(Catch)
catch_discover_tests(controller_sim_test)