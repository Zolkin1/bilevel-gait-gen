
add_executable(demo demo.cpp)

# simulator
if(SIMULATION)
    target_link_libraries(demo simulator)
endif ()

# yaml-cpp
find_package(yaml-cpp REQUIRED)
target_link_libraries(demo yaml-cpp::yaml-cpp)
