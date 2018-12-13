set(ENTRYPOINT_TEST ${TEST_DIRECTORY}/main.cpp)
add_executable(symmetry_tests ${ENTRYPOINT_TEST} ${SRC_DIRECTORY}/signature.cpp
        ${SRC_DIRECTORY}/transformation.cpp)
target_link_libraries(symmetry_tests igl::core igl::opengl_glfw ms-clustering flann_cpp)
