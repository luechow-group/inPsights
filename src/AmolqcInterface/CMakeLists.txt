cmake_minimum_required(VERSION 3.10)

set(NORULES YES)
set(PFUNIT OFF)
add_subdirectory(Amolqc)

add_subdirectory(AmolqcElocInterface)
add_subdirectory(AmolqcFileImport)
