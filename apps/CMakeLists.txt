
add_executable(demo demo.cpp)

# simulator
if(SIMULATION)
    target_link_libraries(demo simulator)
endif ()

# mpc
target_link_libraries(demo mpc)
