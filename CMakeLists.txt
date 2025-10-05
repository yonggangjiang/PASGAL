cmake_minimum_required(VERSION 3.16)

project(PASGAL LANGUAGES C CXX)

# Options
option(PASGAL_ENABLE_WERROR "Enable -Werror on builds" ON)
option(PASGAL_ENABLE_NATIVE_ARCH "Enable -march=native optimization" ON)
option(PASGAL_ENABLE_PARLAY_SEQUENTIAL "Build sequential variants with PARLAY_SEQUENTIAL" OFF)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Common compile flags from Makefiles: -std=c++20 -Wall -Wextra -Werror
add_library(pasgal_warnings INTERFACE)
target_compile_features(pasgal_warnings INTERFACE cxx_std_20)
target_compile_options(pasgal_warnings INTERFACE
  $<$<COMPILE_LANGUAGE:CXX>:-Wall;-Wextra>
)
if(PASGAL_ENABLE_WERROR)
  target_compile_options(pasgal_warnings INTERFACE $<$<COMPILE_LANGUAGE:CXX>:-Werror>)
endif()

# Suppress specific external header warnings that surface in our TU builds
target_compile_options(pasgal_warnings INTERFACE
  $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>>:-Wno-tautological-compare;-Wno-error=tautological-compare>
)

# Includes
set(PASGAL_ROOT ${CMAKE_CURRENT_SOURCE_DIR})

add_library(pasgal_includes INTERFACE)
target_include_directories(pasgal_includes INTERFACE
  ${PASGAL_ROOT}/src
)

# Single-file headers under external/ (e.g., external/connectivity.h)
add_library(external_headers INTERFACE)
target_include_directories(external_headers INTERFACE SYSTEM
  ${PASGAL_ROOT}/external
)

if(PASGAL_ENABLE_NATIVE_ARCH)
  add_library(pasgal_native INTERFACE)
  target_compile_options(pasgal_native INTERFACE -march=native)
else()
  add_library(pasgal_native INTERFACE)
endif()

add_library(pasgal_common INTERFACE)
target_link_libraries(pasgal_common INTERFACE pasgal_warnings pasgal_includes pasgal_native external_headers)

# Helper to create executables succinctly
function(pasgal_add_executable target)
  add_executable(${target} ${ARGN})
  target_link_libraries(${target} PRIVATE pasgal_common)
  target_compile_definitions(${target} PUBLIC QUILL_ROOT_LOGGER_ONLY)
endfunction()

# ========= Subprojects =========

add_subdirectory(external/parlaylib)
set(QUILL_MASTER_PROJECT OFF CACHE BOOL "" FORCE)
add_subdirectory(external/quill)

# Link third-party libraries globally
target_link_libraries(pasgal_common INTERFACE parlay quill)

add_subdirectory(src)

# Convenience default target (aggregate of all executables)
add_custom_target(all_apps)
# Depend on existing targets only
foreach(tgt IN ITEMS
  dinic push-relabel sssp dijkstra bfs seq-bfs bfs_test scc tarjan reachability
  fast-bcc tarjan-vishkin hopcroft-tarjan basic_analytics
  symmetrize pbbs2bin validator generate_random_graph generate_grid_graph get_diameter)
  if(TARGET ${tgt})
    add_dependencies(all_apps ${tgt})
  endif()
endforeach()

message(STATUS "PASGAL configured. Use target 'all_apps' to build all executables.")


