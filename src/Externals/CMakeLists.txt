cmake_minimum_required(VERSION 3.11)

add_subdirectory(googletest)
add_subdirectory(dualmc)
set(YAML_CPP_BUILD_TESTS OFF CACHE BOOL "Don't build the YAML tests.")
set(YAML_CPP_BUILD_CONTRIB OFF CACHE BOOL "Don't build YAML contributions.")
add_subdirectory(yaml-cpp)
add_subdirectory(spdlog)