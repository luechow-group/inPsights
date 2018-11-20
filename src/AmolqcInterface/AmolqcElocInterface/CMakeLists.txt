cmake_minimum_required(VERSION 3.10)

project(AmolqcElocInterface C CXX Fortran)

set ( PROJECT_FOLDER_NAME ${PROJECT_NAME})
set ( PROJECT_PATH ${PROJECT_NAME})
set ( OUTPUT_NAME "${PROJECT_FOLDER_NAME}" )
message ( "-- ++ CREATING PROJECT OR MAKEFILE ---> ${OUTPUT_NAME}" )

set(HEADER_FILES
        include/ElectronicWaveFunction.h
        include/ElectronicWaveFunctionProblem.h
        )

set(SOURCE_FILES
        source/ElectronicWaveFunction.cpp
        source/ElectronicWaveFunctionProblem.cpp
        )

add_library(${PROJECT_NAME} STATIC ${SOURCE_FILES} ${HEADER_FILES})

target_include_directories(${PROJECT_NAME}
        PUBLIC
            include
        PRIVATE
            source
        )

target_link_libraries(${PROJECT_NAME}
        BaseLib
        cppoptlib
        inPsightsInterface
        AmolqcFileImport
        )

#add_gtest_subdirectory_with_dependency(tests "gtest;pthread")
