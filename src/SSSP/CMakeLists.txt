add_executable(sssp ${CMAKE_CURRENT_LIST_DIR}/sssp.cpp)
target_link_libraries(sssp PRIVATE pasgal_common)

add_executable(dijkstra ${CMAKE_CURRENT_LIST_DIR}/dijkstra.cpp)
target_link_libraries(dijkstra PRIVATE pasgal_common)

