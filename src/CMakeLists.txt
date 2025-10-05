add_subdirectory(Max-Flow)

# Build other algorithms/tools only when explicitly enabled
option(PASGAL_BUILD_EXTRA_APPS "Build extra apps beyond Max-Flow" OFF)
if(PASGAL_BUILD_EXTRA_APPS)
  add_subdirectory(SSSP)
  add_subdirectory(BFS)
  add_subdirectory(SCC)
  add_subdirectory(Reachability)
  add_subdirectory(BCC)
  add_subdirectory(basic_analytics)
  add_subdirectory(utils)
endif()

