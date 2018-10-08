cmake_minimum_required(VERSION 3.5)

add_subdirectory(googletest)
add_subdirectory(BSplines)
add_dependencies(AmolqcppTests BSplines_tests.exe)

add_subdirectory(cppoptlib)
add_dependencies(AmolqcppTests cppoptlib_tests.exe)

option(YAML_CPP_BUILD_TESTS "Don't build the YAML tests." OFF)
option(YAML_CPP_BUILD_CONTRIB "Don't build YAML contributions." OFF)
add_subdirectory(yaml-cpp)

add_subdirectory(spdlog)