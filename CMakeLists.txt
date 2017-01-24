cmake_minimum_required(VERSION 3.6)

set(CMAKE_C_COMPILER "/usr/local/bin/gcc-6")
set(CMAKE_CXX_COMPILER "/usr/local/bin/g++-6")
set(CMAKE_Fortran_COMPILER "/usr/local/bin/gfortran ")
set(CMAKE_CXX_STANDARD 11)

#enable_language (C CXX Fortran)

project(AmolqcGUI)
message("C compiler: ${CMAKE_C_COMPILER}")
message("CXX compiler: ${CMAKE_CXX_COMPILER}")
message("Fortran compiler: ${CMAKE_Fortran_COMPILER}")
# DOES NOT WORK - BOOST NEEDED?
#find_package (MKL)
#option (USE_MKL_LIBRARY "Perform linear algebra operations with the MKL library" ${MKL_FOUND})


add_subdirectory(src)
