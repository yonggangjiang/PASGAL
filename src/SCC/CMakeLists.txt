pasgal_add_executable(scc ${CMAKE_CURRENT_LIST_DIR}/scc.cpp)

pasgal_add_executable(tarjan ${CMAKE_CURRENT_LIST_DIR}/tarjan.cpp)

if(PASGAL_ENABLE_PARLAY_SEQUENTIAL)
  pasgal_add_executable(scc_seq ${CMAKE_CURRENT_LIST_DIR}/scc.cpp)
  target_compile_definitions(scc_seq PRIVATE PARLAY_SEQUENTIAL)
  target_compile_options(scc_seq PRIVATE -O3 -mcx16)
endif()

