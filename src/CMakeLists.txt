set(ENTRYPOINT ${SRC_DIRECTORY}/main.cpp Signature.cpp Signature.h)
add_executable(symmetry ${ENTRYPOINT})
target_link_libraries(symmetry igl::core igl::core igl::opengl_glfw ms-clustering flann_cpp)
