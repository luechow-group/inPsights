cmake_minimum_required(VERSION 3.5)


project(SphericalHarmonics CXX)

set ( PROJECT_FOLDER_NAME ${PROJECT_NAME})

set ( PROJECT_PATH ${PROJECT_NAME})
set ( OUTPUT_NAME "${PROJECT_FOLDER_NAME}" )
message ( "-- ++ CREATING PROJECT OR MAKEFILE ---> ${OUTPUT_NAME}" )

set(CMAKE_CXX_STANDARD 14)
add_definitions(
        -O3
        -Wall
        -pedantic
        -pthread
        )

find_package(Eigen3 REQUIRED)
include_directories(
        ${EIGEN3_INCLUDE_DIR}
        spherical-harmonics/sh
        )

set(HEADER_FILES
        spherical-harmonics/sh/spherical_harmonics.h
        spherical-harmonics/sh/default_image.h
        spherical-harmonics/sh/image.h
        )
set(SOURCE_FILES
        spherical-harmonics/sh/spherical_harmonics.cc
        spherical-harmonics/sh/default_image.cc
        )

add_library(${PROJECT_NAME} STATIC ${HEADER_FILES} ${SOURCE_FILES})

target_include_directories(${PROJECT_NAME} PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/spherical-harmonics>
        )
