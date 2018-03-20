cmake_minimum_required(VERSION 3.5)


project(GaussKronrodCartesianIntegration CXX)

set ( PROJECT_FOLDER_NAME ${PROJECT_NAME})

set ( PROJECT_PATH ${PROJECT_NAME})
set ( OUTPUT_NAME "${PROJECT_FOLDER_NAME}" )
message ( "-- ++ CREATING PROJECT OR MAKEFILE ---> ${OUTPUT_NAME}" )

set(CMAKE_CXX_STANDARD 14)

add_definitions(
        -O3
        -DAPL=0
        -DIBM=0
        -DLIN=1
        -DNDEBUG
        -DCOMPUTE_NI_WEIGHTS_ON_THE_FLY

        -Wall
        -Wcast-align
        -Wchar-subscripts
        -Wconversion
        -Wdouble-promotion
        -Wextra
        -Wfatal-errors
        -Wformat-security
        -Wlogical-op
        -Wno-long-long
        -Wno-psabi
        -Wno-variadic-macros
        #-Wnon-virtual-dtor
        -Wpointer-arith
        #-Wundef
        -Wuninitialized
        -Wunused-local-typedefs
        -Wwrite-strings

        -fstrict-aliasing
        -fno-check-new
        -fno-common
        -fvisibility=hidden
        -fPIC

        -c
        -g
        -march=native
        -pedantic
        -pthread
)

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})


add_library(${PROJECT_NAME} INTERFACE)

target_include_directories(${PROJECT_NAME} INTERFACE /
        MPFR/
        include/
        )

add_subdirectory(tests)
