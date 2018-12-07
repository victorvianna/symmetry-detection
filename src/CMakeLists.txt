set(ENTRYPOINT ${SRC_DIRECTORY}/main.cpp)
add_executable(symmetry ${ENTRYPOINT})
target_link_libraries(symmetry igl::core)
