add_executable(fast-bcc ${CMAKE_CURRENT_LIST_DIR}/fast-bcc.cpp)
target_link_libraries(fast-bcc PRIVATE pasgal_common)

add_executable(tarjan-vishkin ${CMAKE_CURRENT_LIST_DIR}/tarjan-vishkin.cpp)
target_link_libraries(tarjan-vishkin PRIVATE pasgal_common)

add_executable(hopcroft-tarjan ${CMAKE_CURRENT_LIST_DIR}/hopcroft-tarjan.cpp)
target_link_libraries(hopcroft-tarjan PRIVATE pasgal_common)

