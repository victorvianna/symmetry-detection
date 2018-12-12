set(CMAKE_CXX_STANDARD 11)
cmake_minimum_required(VERSION 3.1)

project(symmetry-detection)

set(ROOT_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set(SRC_DIRECTORY ${ROOT_DIRECTORY}/src)
set(TEST_DIRECTORY ${ROOT_DIRECTORY}/test)
set(LIB_DIRECTORY ${ROOT_DIRECTORY}/lib)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${ROOT_DIRECTORY}/cmake)

## libigl
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)
option(LIBIGL_WITH_VIEWER            "Use OpenGL viewer"  ON)
find_package(LIBIGL REQUIRED QUIET)

# Mean shift clustering lib
add_subdirectory(${LIB_DIRECTORY}/ms-clustering)
include_directories(${LIB_DIRECTORY}/ms-clustering)

# Nearest neighbors lib
add_subdirectory(${LIB_DIRECTORY}/flann)
include_directories(${LIB_DIRECTORY}/flann/src/cpp)

include_directories(${SRC_DIRECTORY})

add_subdirectory(${TEST_DIRECTORY})
add_subdirectory(${SRC_DIRECTORY})
