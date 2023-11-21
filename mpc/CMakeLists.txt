
set(MPC_HEADER_LIST "include/mpc.h"
        "include/qp_interface.h"
        "include/centroidal_model.h"
        "include/integrator.h"
        "include/trajectory.h"
        "include/spline.h"
        "include/inputs.h"
        "include/euler_integrator.h"
        "include/qp_data.h"
        "include/osqp_interface.h")

add_library(mpc mpc.cpp
        qp_interface.cpp
        centroidal_model.cpp
        integrator.cpp
        trajectory.cpp
        spline.cpp
        inputs.cpp
        euler_integrator.cpp
        osqp_interface.cpp
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

# OSPQ and OSQPEigen
set(${PROFILING} 1)
find_package(osqp REQUIRED)
find_package(OsqpEigen REQUIRED)
target_link_libraries(mpc PUBLIC osqp::osqp OsqpEigen::OsqpEigen)

# MKL and OpenMP
#find_package(OpenMP)
#target_link_libraries(mpc PUBLIC OpenMP::OpenMP_CXX)
#set(THREADS_PREFER_PTHREAD_FLAG ON)
#find_package(Threads REQUIRED)
##target_link_libraries(my_app PRIVATE Threads::Threads)
#
#find_package(MKL CONFIG REQUIRED)
#message(STATUS "${MKL_IMPORTED_TARGETS}") #Provides available list of targets based on input
#add_executable(myapp app.c)
#
#target_compile_options(mpc PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
#target_include_directories(mpc PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
#target_link_libraries(mpc PUBLIC $<LINK_ONLY:MKL::MKL> OpenMP::OpenMP_CXX Threads::Threads)
