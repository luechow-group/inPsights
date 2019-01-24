cmake_minimum_required(VERSION 3.11)

project(DensityBasedScan CXX)

set ( PROJECT_FOLDER_NAME ${PROJECT_NAME} )
set ( PROJECT_PATH ${PROJECT_NAME})
set ( OUTPUT_NAME "${PROJECT_FOLDER_NAME}" )
message ( "-- ++ CREATING PROJECT OR MAKEFILE ---> ${OUTPUT_NAME}" )

set(HEADER_FILES
        include/DensityBasedScan.h
        include/VantagePointTree.h
        )

set(SOURCE_FILES
        source/Helper.cpp
        )

add_library(${PROJECT_NAME} STATIC ${SOURCE_FILES} ${HEADER_FILES})
target_include_directories(${PROJECT_NAME}
        PUBLIC include
        PRIVATE source
        )

target_link_libraries(${PROJECT_NAME}
        spdlog
        StructuralSimilarity
        )

add_gtest_subdirectory_with_dependency(tests "gtest;pthread")


