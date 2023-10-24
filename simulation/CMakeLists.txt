# Setup Mujoco things

# Setup Mujoco path
set(MUJOCO_VERSION "mujoco-3.0.0" CACHE STRING "mujoco version")
set(MUJOCO_DIR "/home/$ENV{USER}/${MUJOCO_VERSION}" CACHE PATH "Path to Mujoco dir")
message(STATUS "mujoco directory: ${MUJOCO_DIR}")

# Setup the Mujoco library
find_library(MUJOCO_LIB mujoco REQUIRED HINTS ${MUJOCO_DIR}/lib)

add_library(mujoco SHARED IMPORTED GLOBAL)
target_include_directories(mujoco INTERFACE ${MUJOCO_DIR}/include)
target_include_directories(mujoco INTERFACE ${MUJOCO_DIR}/include/mujoco)

set_property(TARGET mujoco PROPERTY IMPORTED_LOCATION ${MUJOCO_DIR}/lib/libmujoco.so)

# Create the simulator library

set(SIM_HEADER_LIST "include/simulator.h"
                    "include/simulation_robot.h"
                    "include/controller.h"
                    "include/glfw_adapter.h"
                    "include/array_safety.h"
                    "include/glfw_corevideo.h"
                    "include/glfw_dispatch.h"
                    "include/platform_ui_adapter.h"
                    "include/simulate.h")

add_library(simulator simulator.cpp
                      simulation_robot.cpp
                      controller.cpp
                    glfw_adapter.cc
                    glfw_dispatch.cc
                    platform_ui_adapter.cc
                    simulate.cc
                      ${SIM_HEADER_LIST})

target_include_directories(simulator PUBLIC include)
target_link_libraries(simulator PUBLIC mujoco)

# GLFW
find_package(glfw3 REQUIRED)
target_link_libraries(simulator PUBLIC glfw)

# Eigen
find_package(Eigen3 REQUIRED)
target_link_libraries(simulator PUBLIC Eigen3::Eigen)

# yaml-cpp
find_package(yaml-cpp REQUIRED)
target_link_libraries(simulator PRIVATE yaml-cpp::yaml-cpp)

