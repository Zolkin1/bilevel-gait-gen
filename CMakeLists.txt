cmake_minimum_required(VERSION 3.16)

project(bilevel-gait-gen
    VERSION 0.1
    LANGUAGES CXX)

# TODO: add find packages

set(SIMULATION ON CACHE BOOL "Configure if the simulator is in use.")

if(SIMULATION)
    add_subdirectory(simulation)
endif ()

add_subdirectory(apps)

