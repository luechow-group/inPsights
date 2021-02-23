cmake_minimum_required(VERSION 3.11)

option(BUILD_GUI    "ON/OFF - Do/Don't build the graphical user interface." OFF)
option(BUILD_EIGEN  "ON/OFF - Build/Try to find the required Eigen3 library." OFF)
option(YAML_CPP_BUILD_TESTS "ON/OFF - yaml test" OFF)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake/Modules")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(inPsights)

include(cmake/findSharedMemoryParallelismLibraries.cmake)
include(cmake/buildInfo.cmake)

if(DEFINED ENV{BOOST_ROOT})
    set(BOOST_ROOT "$ENV{BOOST_ROOT}")
    message("-- BOOST_ROOT is set: ${BOOST_ROOT}")
else(DEFINED ENV{BOOST_ROOT})
    message("-- Note: BOOST_ROOT not set. Let's see whether we can find BOOST anyway.")
endif(DEFINED ENV{BOOST_ROOT})
find_package(Boost)
include_directories(${Boost_INCLUDE_DIRS})

if (${BUILD_EIGEN})
    message("BUILD_EIGEN=${BUILD_EIGEN} Building Eigen3...")
    include(cmake/downloadEigen.cmake)
else()
    message("BUILD_EIGEN=${BUILD_EIGEN} Trying to find Eigen3...")
    find_package(Eigen3 REQUIRED)
    include_directories(${EIGEN3_INCLUDE_DIR})
endif()

include(cmake/downloadGoogletest.cmake)
include(cmake/downloadSpdlog.cmake)
include(cmake/downloadYamlCpp.cmake)


# putting all executables in cmake-build-*/bin
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/bin)

# Ctest related
enable_testing()
add_custom_target(inPsightsTests COMMAND "")
include(cmake/addGtestSubdirectoryWithDependency.cmake)
add_subdirectory(src)
