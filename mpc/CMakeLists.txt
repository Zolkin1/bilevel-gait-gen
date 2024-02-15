
set(MPC_HEADER_LIST "include/mpc.h"
        "include/qp/qp_interface.h"
        "include/models/centroidal_model.h"
        "include/integrator.h"
        "include/trajectory.h"
        "include/spline/spline.h"
        "include/euler_integrator.h"
        "include/qp/qp_data.h"
        "include/qp/osqp_interface.h"
        "include/rk_integrator.h"
        "include/gait_optimizer.h"
        "include/models/single_rigid_body_model.h"
        "include/models/model.h"
        "include/mpc_single_rigid_body.h"
        "include/mpc_centroidal.h"
        "include/spline/spline_node.h"
        "include/spline/end_effector_splines.h"
        "include/qp/clarabel_interface.h")

add_library(mpc mpc.cpp
        qp/qp_interface.cpp
        models/centroidal_model.cpp
        integrator.cpp
        trajectory.cpp
        spline/spline.cpp
        euler_integrator.cpp
        qp/osqp_interface.cpp
        rk_integrator.cpp
        qp/qp_data.cpp
        gait_optimizer.cpp
        models/single_rigid_body_model.cpp
        models/model.cpp
        mpc_single_rigid_body.cpp
        mpc_centroidal.cpp
        spline/spline_node.cpp
        spline/end_effector_splines.cpp
        qp/clarabel_interface.cpp
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

# Clarabel
add_library(clarabel STATIC IMPORTED)
set_property(TARGET clarabel PROPERTY
        IMPORTED_LOCATION "/home/zolkin/Repos/Clarabel.cpp/rust_wrapper/target/release/libclarabel_c.so")
target_include_directories(mpc PUBLIC "/home/zolkin/Repos/Clarabel.cpp/include")
target_link_libraries(mpc PUBLIC clarabel)

#lsqr
message(STATUS "CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}")
find_package(lsqr-cpp REQUIRED)
target_link_libraries(mpc PUBLIC lsqr-cpp::lsqr-cpp)
