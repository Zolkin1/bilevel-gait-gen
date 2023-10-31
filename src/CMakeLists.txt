# Source cmake

set(MPC_HEADER_LIST "../include/config_parser.h")

add_library(mpc config_parser.cpp ${MPC_HEADER_LIST})
target_include_directories(mpc PUBLIC ../include)

# Eigen
find_package(Eigen3 REQUIRED)
target_link_libraries(mpc PUBLIC Eigen3::Eigen)

# yaml-cpp
find_package(yaml-cpp REQUIRED)
target_link_libraries(mpc PRIVATE yaml-cpp::yaml-cpp)