cmake_minimum_required(VERSION 3.5)

if (UNIX)
    if(APPLE)
        message("using apple compiler options")
        set(CMAKE_C_COMPILER "/usr/local/bin/gcc-7")
        set(CMAKE_CXX_COMPILER "/usr/local/bin/g++-7")
        set(CMAKE_Fortran_COMPILER "/usr/local/bin/gfortran-7")
    else()
        message("using UNIX compiler options")
        set(CMAKE_C_COMPILER gcc)
        set(CMAKE_CXX_COMPILER g++)
        set(CMAKE_Fortran_COMPILER gfortran)
    endif(APPLE)
endif(UNIX)

set(CMAKE_CXX_STANDARD 11)

project(AmolqcGUI)

message("C compiler: ${CMAKE_C_COMPILER}")
message("CXX compiler: ${CMAKE_CXX_COMPILER}")
message("Fortran compiler: ${CMAKE_Fortran_COMPILER}")

find_package(OpenMP REQUIRED)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")

#find_package (MKL)
#option (USE_MKL_LIBRARY "Perform linear algebra operations with the MKL library" ${MKL_FOUND})

# Prevent error from policy CMP0042. See >>> cmake  --help-policy CMP0042
set(CMAKE_MACOSX_RPATH 1)

#------------------------------------------------#
# Download and unpack googletest at configure time
configure_file(DownloadGoogleTest.in
        googletest-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download )
execute_process(COMMAND ${CMAKE_COMMAND} --build .
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download )

# Prevent GoogleTest from overriding our compiler/linker options
# when building with Visual Studio
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Add googletest directly to our build. This adds
# the following targets: gtest, gtest_main, gmock
# and gmock_main

add_subdirectory(${CMAKE_BINARY_DIR}/googletest-src
        ${CMAKE_BINARY_DIR}/googletest-build)

# The gtest/gmock targets carry header search path
# dependencies automatically when using CMake 2.8.11 or
# later. Otherwise we have to add them here ourselves.
if (CMAKE_VERSION VERSION_LESS 2.8.11)
    include_directories(
            "${gtest_SOURCE_DIR}/include"
            "${gmock_SOURCE_DIR}/include")
endif()


##---------------------------------------------------#
## Download and unpack Eigen Library at configure time
#configure_file(DownloadEigen.in
#        eigen-download/CMakeLists.txt)
#execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
#        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/eigen-download )
#execute_process(COMMAND ${CMAKE_COMMAND} --build .
#        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/eigen-download )
#
## Add Eigen library directly to our build. This adds
## the following targets: eigen
#
#add_subdirectory(
#        ${CMAKE_BINARY_DIR}/eigen-src
#        ${CMAKE_BINARY_DIR}/eigen-build
#        )
find_package(Eigen3 REQUIRED)
include_directories( ${EIGEN3_INCLUDE_DIR} )

add_subdirectory(src)
