cmake_minimum_required(VERSION 3.5)

include(cmake/SetCompilers.cmake)

set(CMAKE_CXX_STANDARD 14)

project(Amolqcpp)

# Prevent error from policy CMP0042. See >>> cmake  --help-policy CMP0042
set(CMAKE_MACOSX_RPATH 1)

include(cmake/FindSharedMemoryParallelismLibraries.cmake)

include(cmake/DownloadGoogleTest.cmake)

# switch between external eigen and cloning it from github
#include(cmake/DownloadEigen.cmake)
find_package(Eigen3 REQUIRED)
include_directories( ${EIGEN3_INCLUDE_DIR} )

add_subdirectory(src)
