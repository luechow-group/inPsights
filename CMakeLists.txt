cmake_minimum_required(VERSION 3.5)
project(LocalSpinMultiplicity)

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR})
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

include_directories(include)

set(SOURCE_FILES
        src/main.cpp
        include/pse.h
        src/Particle.cpp
        include/Particle.h
        src/Electron.cpp
        src/Electron.h
        src/Core.cpp
        src/Core.h
    )

add_executable(LocalSpinMultiplicity ${SOURCE_FILES})
