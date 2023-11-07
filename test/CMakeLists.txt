find_package(Catch2 3 REQUIRED)

add_executable(sim_test simulator_test.cpp)

target_link_libraries(sim_test PRIVATE Catch2::Catch2WithMain)

target_link_libraries(sim_test PUBLIC simulator)
target_link_libraries(sim_test PUBLIC controller)
target_link_libraries(sim_test PUBLIC utils)

include(Catch)
catch_discover_tests(sim_test)