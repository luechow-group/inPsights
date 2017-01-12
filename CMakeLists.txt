cmake_minimum_required(VERSION 3.6)

project(AmolqcGUI)
message("C compiler: ${CMAKE_C_COMPILER}")
message("CXX compiler: ${CMAKE_CXX_COMPILER}")
message("Fortran compiler: ${CMAKE_Fortran_COMPILER}")
# DOES NOT WORK - BOOST NEEDED?
#find_package (MKL)
#option (USE_MKL_LIBRARY "Perform linear algebra operations with the MKL library" ${MKL_FOUND})


add_subdirectory(src)
