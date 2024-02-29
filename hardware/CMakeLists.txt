#set(PINOCCHIO_DIR "/home/$ENV{USER}/repos/pinocchio" CACHE PATH "Path to pinocchio")
#set(CMAKE_MODULE_PATH "${PINOCCHIO_DIR}/cmake/find-external/CppAD")

# Create the controller library
set(HARDWARE_HEADER_LIST
        "include/hardware_robot.h"
        "unitree_lib/udp.h")

add_library(hardware
        hardware_robot.cpp
        ${HARDWARE_HEADER_LIST})

add_executable(hardware_interface hardware_interface.cpp)
target_link_libraries(hardware_interface PUBLIC hardware)
target_link_libraries(hardware PUBLIC /home/zolkin/AmberLab/bilevel-gait-gen/hardware/unitree_lib/libunitree_legged_sdk_amd64.so)

target_include_directories(hardware PUBLIC include unitree_lib)

add_executable(optitrack_client optitrack_client.cpp)

# Eigen
find_package(Eigen3 REQUIRED)
target_link_libraries(hardware PUBLIC Eigen3::Eigen)

# MPC
target_link_libraries(hardware PUBLIC mpc)

# Controller
target_link_libraries(hardware PUBLIC controller)

# Utils
target_link_libraries(hardware_interface PUBLIC utils)

# LCM
target_link_libraries(hardware PUBLIC lcm)

# optitrack stream client
find_package(optitrack-stream-client REQUIRED)
target_link_libraries(hardware PUBLIC optitrack-stream-client::stream_client)
target_link_libraries(optitrack_client PUBLIC optitrack-stream-client::stream_client)