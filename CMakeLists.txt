cmake_minimum_required(VERSION 3.5)

include(cmake/SetCompilers.cmake)

set(CMAKE_CXX_STANDARD 14)

project(Amolqcpp)

# Prevent error from policy CMP0042. See >>> cmake  --help-policy CMP0042
set(CMAKE_MACOSX_RPATH 1)

include(cmake/FindSharedMemoryParallelismLibraries.cmake)

include(cmake/DownloadGoogleTest.cmake)

if (DEFINED BUILD_EIGEN AND (${BUILD_EIGEN}))
    message("BUILD_EIGEN build eigen is ON. Building Eigen3...")
    include(cmake/DownloadEigen.cmake)
else()
    message("BUILD_EIGEN build eigen is undefined or OFF. Trying to find Eigen3...")
    find_package(Eigen3 REQUIRED)
    include_directories(${EIGEN3_INCLUDE_DIR})
endif()

add_subdirectory(src)
