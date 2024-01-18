
set(MPC_HEADER_LIST "include/mpc.h"
        "include/qp_interface.h"
        "include/centroidal_model.h"
        "include/integrator.h"
        "include/trajectory.h"
        "include/spline.h"
        "include/inputs.h"
        "include/euler_integrator.h"
        "include/qp_data.h"
        "include/osqp_interface.h"
        "include/rk_integrator.h"
        "gait_optimizer.cpp")

add_library(mpc mpc.cpp
        qp_interface.cpp
        centroidal_model.cpp
        integrator.cpp
        trajectory.cpp
        spline.cpp
        inputs.cpp
        euler_integrator.cpp
        osqp_interface.cpp
        rk_integrator.cpp
        qp_data.cpp
        gait_optimizer.cpp
        ${MPC_HEADER_LIST})
target_include_directories(mpc PUBLIC include)

# Eigen
find_package(Eigen3 REQUIRED)
target_link_libraries(mpc PUBLIC Eigen3::Eigen)

# pinocchio
# TODO: make not hard coded
set(PINOCCHIO_DIR "/home/$ENV{USER}/repos/pinocchio" CACHE PATH "Path to pinocchio")
set(CMAKE_MODULE_PATH "${PINOCCHIO_DIR}/cmake/find-external/CppAD")
find_package(pinocchio)
target_link_libraries(mpc PRIVATE pinocchio)


# TODO: Can go back and play with MKL more, QLDL seems quicker by a noticable amount
# MKL (for OSQP)
set(MKL_DIR /opt/intel/oneapi/mkl/latest/lib/cmake/mkl)
set(MKL_THREADING gnu_thread)
set(MKL_ARCH intel64)
set(MKL_MPI openmpi)
find_package(MKL CONFIG REQUIRED)

# OpenMP (for MKL)
find_package(OpenMP REQUIRED)
target_link_libraries(mpc PUBLIC OpenMP::OpenMP_CXX)

# OSPQ and OSQPEigen
find_package(osqp REQUIRED)
find_package(OsqpEigen REQUIRED)
target_link_libraries(mpc PUBLIC osqp::osqpstatic OsqpEigen::OsqpEigen)

# Controller
target_link_libraries(mpc PUBLIC controller)

# Utils
target_link_libraries(mpc PUBLIC utils)
