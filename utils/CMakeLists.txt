
set(UTIL_HEADER_LIST "include/config_parser.h"
        "include/logger.h"
        "include/timer.h")

add_library(utils config_parser.cpp
        logger.cpp
        timer.cpp
        ${UTIL_HEADER_LIST})
target_include_directories(utils PUBLIC include)

# Eigen
find_package(Eigen3 REQUIRED)
target_link_libraries(utils PUBLIC Eigen3::Eigen)

# yaml-cpp
find_package(yaml-cpp REQUIRED)
target_link_libraries(utils PRIVATE yaml-cpp::yaml-cpp)