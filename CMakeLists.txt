cmake_minimum_required(VERSION 3.16)

project(bilevel-gait-gen
    VERSION 0.1
    LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17 CACHE STRING "The C++ standard to use")


# TODO: add find packages

set(SIMULATION ON CACHE BOOL "Configure if the simulator is in use.")

if (SIMULATION)
    add_subdirectory(simulation)
endif ()

if (CMAKE_PROJECT_NAME STREQUAL bilevel-gait-gen)
    include(CTest)
endif ()

if(CMAKE_PROJECT_NAME STREQUAL bilevel-gait-gen AND BUILD_TESTING)
    add_subdirectory(test)
endif()

add_subdirectory(utils)

add_subdirectory(mpc)

add_subdirectory(controllers)

add_subdirectory(apps)

