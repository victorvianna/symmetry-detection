set(ENTRYPOINT ${SRC_DIRECTORY}/main.cpp signature.cpp)
add_executable(symmetry ${ENTRYPOINT})
target_link_libraries(symmetry igl::core igl::core igl::opengl_glfw ms-clustering flann_cpp)
