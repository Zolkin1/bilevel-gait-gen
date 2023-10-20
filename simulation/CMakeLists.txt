# Setup Mujoco things

# Setup Mujoco path
set(MUJOCO_VERSION "mujoco-3.0.0" CACHE STRING "mujoco version")
set(MUJOCO_DIR "/home/$ENV{USER}/${MUJOCO_VERSION}" CACHE PATH "Path to Mujoco dir")
message(STATUS "mujoco directory: ${MUJOCO_DIR}")

# Setup the Mujoco library
find_library(MUJOCO_LIB mujoco REQUIRED HINTS ${MUJOCO_DIR}/lib)

add_library(mujoco SHARED IMPORTED GLOBAL)
target_include_directories(mujoco INTERFACE ${MUJOCO_DIR}/include)
target_include_directories(mujoco INTERFACE ${MUJOCO_DIR}/include/mujoco)

set_property(TARGET mujoco PROPERTY IMPORTED_LOCATION ${MUJOCO_DIR}/lib/libmujoco.so)

# Create the simulator library
set(SIM_HEADER_LIST "include/simulator.h")

add_library(simulator simulator.cpp ${SIM_HEADER_LIST})
target_include_directories(simulator PUBLIC include)
target_link_libraries(simulator PUBLIC mujoco)