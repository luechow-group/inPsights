cmake_minimum_required(VERSION 3.6)


#============= PROJECT LibAmolqc ===================
project(LibAmolqc)

enable_language (C CXX Fortran)

set ( PROJECT_FOLDER_NAME ${PROJECT_NAME})

set ( PROJECT_PATH ${PROJECT_NAME})
set ( OUTPUT_NAME "${PROJECT_FOLDER_NAME}" )
message ( "-- ++ CREATING PROJECT OR MAKEFILE ---> ${OUTPUT_NAME}" )



set(CMAKE_C_COMPILER "/usr/local/bin/gcc-6")

set(CMAKE_CXX_COMPILER "/usr/local/bin/g++-6")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-g -Wall -fopenmp -llapack -lblas -lgfortran")

set(CMAKE_Fortran_COMPILER "/usr/local/bin/gfortran ")
set(CMAKE_Fortran_FLAGS "-g -fopenmp -llapack -lblas -lstdc++")


#set(AMOLQC /Users/michaelheuer/amolqcGUI/src/LibAmolqc/amolqc)
#message("${AMOLQC}")

find_package(Eigen3 REQUIRED)

include_directories(
        include
        source
        ${GTEST_INCLUDE_DIRS}
        ${EIGEN3_INCLUDE_DIR}
        )

set(HEADER_FILES
        include/ElectronicWaveFunction.h
        )

set(SOURCE_FILES
        source/ElectronicWaveFunction.cpp
        )

link_directories(
        amolqc/src
        amolqc/utils
        ) #searches the specified directories for libraries

add_executable("${PROJECT_NAME}_main" ${SOURCE_FILES} ${HEADER_FILES} main.cpp)
target_link_libraries("${PROJECT_NAME}_main"
        amolqc
        utils
        ) #adds links to amolqc/utils/libamolqc.a and amolqc/src/libutils.a

#add_library(${PROJECT_NAME} STATIC ${SOURCE_FILES} ${HEADER_FILES})
#target_include_directories(${PROJECT_NAME} PUBLIC
#        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
#        )
#add_library(AmolqcGUI::${PROJECT_NAME} ALIAS ${PROJECT_NAME})

#add_subdirectory(tests)
