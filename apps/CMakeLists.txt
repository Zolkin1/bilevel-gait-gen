
add_executable(demo demo.cpp)

# simulator
if(SIMULATION)
    target_link_libraries(demo simulator)
endif ()

# controller
target_link_libraries(demo controller)

# utils
target_link_libraries(demo utils)

# mpc
#target_link_libraries(demo mpc)
