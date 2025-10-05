add_executable(symmetrize ${CMAKE_CURRENT_LIST_DIR}/symmetrize.cpp)
target_link_libraries(symmetrize PRIVATE pasgal_common)

add_executable(pbbs2bin ${CMAKE_CURRENT_LIST_DIR}/pbbs2bin.cpp)
target_link_libraries(pbbs2bin PRIVATE pasgal_common)

add_executable(validator ${CMAKE_CURRENT_LIST_DIR}/validator.cpp)
target_link_libraries(validator PRIVATE pasgal_common)

add_executable(generate_random_graph ${CMAKE_CURRENT_LIST_DIR}/generate_random_graph.cpp)
target_link_libraries(generate_random_graph PRIVATE pasgal_common)

add_executable(generate_grid_graph ${CMAKE_CURRENT_LIST_DIR}/generate_grid_graph.cpp)
target_link_libraries(generate_grid_graph PRIVATE pasgal_common)

add_executable(get_diameter ${CMAKE_CURRENT_LIST_DIR}/get_diameter.cpp)
target_link_libraries(get_diameter PRIVATE pasgal_common)

