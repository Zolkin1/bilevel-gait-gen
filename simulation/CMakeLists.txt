# Setup Mujoco things

# Setup Mujoco path
set(MUJOCO_VERSION "mujoco-3.0.0" CACHE STRING "mujoco version")
set(MUJOCO_DIR "/home/$ENV{USER}/${MUJOCO_VERSION}" CACHE PATH "Path to Mujoco dir")
message(STATUS "mujoco directory: ${MUJOCO_DIR}")

# TODO: Fix this
set(PINOCCHIO_DIR "/home/$ENV{USER}/repos/pinocchio" CACHE PATH "Path to pinocchio")
set(CMAKE_MODULE_PATH "${PINOCCHIO_DIR}/cmake/find-external/CppAD")

# Setup the Mujoco library
find_library(MUJOCO_LIB mujoco REQUIRED HINTS ${MUJOCO_DIR}/lib)

add_library(mujoco SHARED IMPORTED GLOBAL)
target_include_directories(mujoco INTERFACE ${MUJOCO_DIR}/include)
target_include_directories(mujoco INTERFACE ${MUJOCO_DIR}/include/mujoco)

set_property(TARGET mujoco PROPERTY IMPORTED_LOCATION ${MUJOCO_DIR}/lib/libmujoco.so)

# Create the simulator library

set(SIM_HEADER_LIST "include/simulator.h"
                    "include/simulation_robot.h"
                    "include/pd_grav_comp.h"
                    "include/controller.h"
                    "include/glfw_adapter.h"
                    "include/array_safety.h"
                    "include/glfw_corevideo.h"
                    "include/glfw_dispatch.h"
                    "include/platform_ui_adapter.h"
                    "include/simulate.h"
                    "include/qp_control.h")

add_library(simulator mujoco_sim/simulator.cpp
                      simulation_robot.cpp
                      controller.cpp
                    pd_grav_comp.cpp
        mujoco_sim/glfw_adapter.cc
        mujoco_sim/glfw_dispatch.cc
        mujoco_sim/platform_ui_adapter.cc
        mujoco_sim/simulate.cc
        qp_control.cpp
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

# pinocchio
find_package(pinocchio)
target_link_libraries(simulator PRIVATE pinocchio)

# OSPQ and OSQPEigen
find_package(osqp REQUIRED)
find_package(OsqpEigen REQUIRED)
target_link_libraries(simulator PUBLIC osqp::osqp OsqpEigen::OsqpEigen)

