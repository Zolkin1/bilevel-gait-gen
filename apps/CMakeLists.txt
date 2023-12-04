
add_executable(demo demo.cpp
        mpc_sim_demo.cpp)
add_executable(mpc_demo mpc_demo.cpp)
add_executable(mpc_sim mpc_sim_demo.cpp)


# simulator
if(SIMULATION)
    target_link_libraries(demo simulator)
    target_link_libraries(mpc_demo simulator)
    target_link_libraries(mpc_sim simulator)
endif ()

# controller
target_link_libraries(demo controller)
target_link_libraries(mpc_demo controller)
target_link_libraries(mpc_sim controller)

# utils
target_link_libraries(demo utils)
target_link_libraries(mpc_demo utils)
target_link_libraries(mpc_sim utils)

# mpc
target_link_libraries(demo mpc)
target_link_libraries(mpc_demo mpc)
target_link_libraries(mpc_sim mpc)

add_custom_command(TARGET mpc_demo POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${CMAKE_CURRENT_SOURCE_DIR}/a1_configuration.yaml
        ${CMAKE_CURRENT_SOURCE_DIR}/adam_configuration.yaml
        $<TARGET_FILE_DIR:mpc_demo>
)
