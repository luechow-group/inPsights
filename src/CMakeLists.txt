
# get the add_to_source_group function
#include(cmake/MyCmakeFunctions.cmake)

project(AmolqcGUI)

enable_language (C CXX Fortran)
set(CMAKE_CXX_STANDARD 11)

#add_subdirectory(amolqc)

add_subdirectory(BSplines)
add_subdirectory(GUI)
add_subdirectory(Analysis)
add_subdirectory(Methods)