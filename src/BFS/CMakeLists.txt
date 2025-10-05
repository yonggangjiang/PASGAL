add_executable(bfs ${CMAKE_CURRENT_LIST_DIR}/bfs.cpp)
target_link_libraries(bfs PRIVATE pasgal_common)

add_executable(seq-bfs ${CMAKE_CURRENT_LIST_DIR}/seq-bfs.cpp)
target_link_libraries(seq-bfs PRIVATE pasgal_common)

add_executable(bfs_test ${CMAKE_CURRENT_LIST_DIR}/bfs_test.cpp)
target_link_libraries(bfs_test PRIVATE pasgal_common)

