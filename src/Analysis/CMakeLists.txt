cmake_minimum_required(VERSION 3.6)


#============= PROJECT Analysis ===================
project(Analysis)

set ( PROJECT_FOLDER_NAME ${PROJECT_NAME} )

set ( PROJECT_PATH ${PROJECT_NAME})
set ( OUTPUT_NAME "${PROJECT_FOLDER_NAME}" )
message ( "-- ++ CREATING PROJECT OR MAKEFILE ---> ${OUTPUT_NAME}" )

set(CMAKE_CXX_STANDARD 11)
#set(CMAKE_CXX_FLAGS "-g -Wall") # TODO needed?


# Find Eigen lib
find_package(Eigen3 REQUIRED)

include_directories(
        include
        source
        ${GTEST_INCLUDE_DIRS}
        ${EIGEN3_INCLUDE_DIR}
        )

set(HEADER_FILES
        include/IWaveFunction.h
        )

set(SOURCE_FILES
        source/IWaveFunction.cpp
        )

add_library(${PROJECT_NAME} STATIC ${SOURCE_FILES} ${HEADER_FILES})


target_include_directories(${PROJECT_NAME} PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        )

add_library(AmolqcGUI::${PROJECT_NAME} ALIAS ${PROJECT_NAME})

add_subdirectory(tests)
