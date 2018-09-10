cmake_minimum_required(VERSION 3.5)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake/Modules")

include(cmake/SetCompilers.cmake)

set(CMAKE_CXX_STANDARD 17)

project(Amolqcpp)

# Prevent error from policy CMP0042. See >>> cmake  --help-policy CMP0042
set(CMAKE_MACOSX_RPATH 1)

message("${CMAKE_MODULE_PATH}")

include(cmake/FindSharedMemoryParallelismLibraries.cmake)

include(cmake/DownloadGoogleTest.cmake)

if(DEFINED ENV{BOOST_ROOT})
    set(BOOST_ROOT "$ENV{BOOST_ROOT}")
    message("-- BOOST_ROOT is set: ${BOOST_ROOT}")
else(DEFINED ENV{BOOST_ROOT})
    message("-- Note: BOOST_ROOT not set. Let's see whether we can find BOOST anyway.")
endif(DEFINED ENV{BOOST_ROOT})
find_package(Boost)
include_directories(${Boost_INCLUDE_DIRS})

if (DEFINED BUILD_EIGEN AND (${BUILD_EIGEN}))
    message("BUILD_EIGEN build eigen is ON. Building Eigen3...")
    include(cmake/DownloadEigen.cmake)
else()
    message("BUILD_EIGEN build eigen is undefined or OFF. Trying to find Eigen3...")
    find_package(Eigen3 REQUIRED)
    include_directories(${EIGEN3_INCLUDE_DIR})
endif()

add_subdirectory(src)
