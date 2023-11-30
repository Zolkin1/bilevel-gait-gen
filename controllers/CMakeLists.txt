# TODO: Fix this (make it not hard coded)
set(PINOCCHIO_DIR "/home/$ENV{USER}/repos/pinocchio" CACHE PATH "Path to pinocchio")
set(CMAKE_MODULE_PATH "${PINOCCHIO_DIR}/cmake/find-external/CppAD")

# Create the controller library
set(CONTROLLER_HEADER_LIST "include/pd_grav_comp.h"
        "include/controller.h"
        "include/qp_control.h"
        "include/mpc_controller.h")

add_library(controller controller.cpp
        pd_grav_comp.cpp
        qp_control.cpp
        mpc_controller.cpp
        ${CONTROLLER_HEADER_LIST})

target_include_directories(controller PUBLIC include)

# Eigen
find_package(Eigen3 REQUIRED)
target_link_libraries(controller PUBLIC Eigen3::Eigen)

# pinocchio
find_package(pinocchio)
target_link_libraries(controller PRIVATE pinocchio)

# OSPQ and OSQPEigen
find_package(osqp REQUIRED)
find_package(OsqpEigen REQUIRED)
target_link_libraries(controller PUBLIC osqp::osqp OsqpEigen::OsqpEigen)

# MPC
target_link_libraries(controller PUBLIC mpc)