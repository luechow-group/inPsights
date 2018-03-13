cmake_minimum_required(VERSION 3.5)

project(LocalPotentialEnergy C CXX Fortran)


set ( PROJECT_FOLDER_NAME ${PROJECT_NAME} )
set ( PROJECT_PATH ${PROJECT_NAME})
set ( OUTPUT_NAME "${PROJECT_FOLDER_NAME}" )
message ( "-- ++ CREATING PROJECT OR MAKEFILE ---> ${OUTPUT_NAME}" )

set(CMAKE_CXX_STANDARD 11)

find_package(Eigen3 REQUIRED)
include_directories( ${EIGEN3_INCLUDE_DIR} )

include_directories(
        include
        source
        resources
)

set(HEADER_FILES
        include/PotentialProblem.h
        include/LagrangeProblem.h
        include/TestProblems.h
        )

set(SOURCE_FILES
        source/PotentialProblem.cpp
        )

file(GLOB_RECURSE RESOURCE_FILES resources/*.*)
file(COPY ${RESOURCE_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_library(${PROJECT_NAME} STATIC ${SOURCE_FILES} ${HEADER_FILES})
target_include_directories(${PROJECT_NAME} PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        )
target_link_libraries(${PROJECT_NAME}
        BaseLib
        AmolqcInterface
        )

add_subdirectory(tests)
