
# get the add_to_source_group function
#include(cmake/MyCmakeFunctions.cmake)

project(AmolqcGUI)

set(CMAKE_CXX_STANDARD 11)

add_subdirectory(BSplines)
add_subdirectory(GUI)
