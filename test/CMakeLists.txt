set(ENTRYPOINT_TEST ${TEST_DIRECTORY}/main.cpp)
add_executable(symmetry_tests ${ENTRYPOINT_TEST})
target_link_libraries(symmetry_tests igl::core igl::opengl_glfw)
