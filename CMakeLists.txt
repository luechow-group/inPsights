cmake_minimum_required(VERSION 3.6)

if (UNIX)
    message("using UNIX compiler options")
    set(CMAKE_C_COMPILER gcc)
    set(CMAKE_CXX_COMPILER g++)
    set(CMAKE_Fortran_COMPILER gfortran)
endif(UNIX)

if(APPLE)
    message("using apple compiler options")
    set(CMAKE_C_COMPILER "/usr/local/bin/gcc-6")
    set(CMAKE_CXX_COMPILER "/usr/local/bin/g++-6")
    set(CMAKE_Fortran_COMPILER "/usr/local/bin/gfortran ")
endif(APPLE)

set(CMAKE_CXX_FLAGS "-g -Wall -fopenmp -llapack -lblas -lgfortran")
set(CMAKE_Fortran_FLAGS "-g -fopenmp -llapack -lblas -lstdc++")
set(CMAKE_CXX_STANDARD 11)

project(AmolqcGUI)
message("C compiler: ${CMAKE_C_COMPILER}")
message("CXX compiler: ${CMAKE_CXX_COMPILER}")
message("Fortran compiler: ${CMAKE_Fortran_COMPILER}")

#find_package (MKL)
#option (USE_MKL_LIBRARY "Perform linear algebra operations with the MKL library" ${MKL_FOUND})


add_subdirectory(src)
