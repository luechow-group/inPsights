cmake_minimum_required(VERSION 3.11)

option(PFUNIT "" OFF)
option(CONFIG_FILE "" OFF)
option(WARNINGS "" OFF)

add_subdirectory(Amolqc)

add_subdirectory(AmolqcElocInterface)
add_subdirectory(AmolqcFileImport)
