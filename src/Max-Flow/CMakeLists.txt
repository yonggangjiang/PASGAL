add_executable(dinic ${CMAKE_CURRENT_LIST_DIR}/dinic.cpp)
target_link_libraries(dinic PRIVATE pasgal_common)
target_compile_definitions(dinic PRIVATE QUILL_ROOT_LOGGER_ONLY)

add_executable(push-relabel ${CMAKE_CURRENT_LIST_DIR}/push-relabel.cpp)
target_link_libraries(push-relabel PRIVATE pasgal_common)
target_compile_definitions(push-relabel PRIVATE QUILL_ROOT_LOGGER_ONLY)

