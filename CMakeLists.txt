cmake_minimum_required(VERSION 3.16)

project(bilevel-gait-gen
    VERSION 0.1
    LANGUAGES CXX)

message(STATUS "CMAKE Build Type: ${CMAKE_BUILD_TYPE}")
message(STATUS "CMAKE_CXX_FLAGS_RELWITHDEBINFO: " ${CMAKE_CXX_FLAGS_RELWITHDEBINFO})
message(STATUS "CMAKE_CXX_FLAGS_DEBUG: " ${CMAKE_CXX_FLAGS_DEBUG})
message(STATUS "CMAKE_CXX_FLAGS_RELEASE: " ${CMAKE_CXX_FLAGS_RELEASE})

set(CMAKE_CXX_STANDARD 17 CACHE STRING "The C++ standard to use")


# TODO: add find packages

#add_compile_options(-ffast-math)

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

add_subdirectory(hardware)

