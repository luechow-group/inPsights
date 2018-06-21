cmake_minimum_required(VERSION 3.5)

add_subdirectory(BSplines)
add_subdirectory(cppoptlib)
add_subdirectory(json)

set(YAML_CPP_BUILD_TESTS OFF)
add_subdirectory(yaml-cpp)
