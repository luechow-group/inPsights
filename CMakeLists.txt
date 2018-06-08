cmake_minimum_required(VERSION 3.5)

include(cmake/SetCompilers.cmake)

set(CMAKE_CXX_STANDARD 14)

project(Amolqcpp)

# Prevent error from policy CMP0042. See >>> cmake  --help-policy CMP0042
set(CMAKE_MACOSX_RPATH 1)

include(cmake/DownloadGoogleTest.cmake)
set(BUILT_GOOGLETEST_ALREADY TRUE)


find_package(OpenMP REQUIRED)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")

find_package(Eigen3 REQUIRED)
include_directories( ${EIGEN3_INCLUDE_DIR} )

add_subdirectory(src)
