cmake_minimum_required(VERSION 3.0.2)
enable_testing()

find_package(GTest REQUIRED)

include_directories(SYSTEM
        ${GTEST_INCLUDE_DIRS}
        )

FILE(GLOB TEST_SOURCE_FILES unit_tests/*)

add_executable("${PROJECT_NAME}_RunAllTests" ${TEST_SOURCE_FILES})

link_directories(
        amolqc/src
        amolqc/utils
) #searches the specified directories for libraries

target_link_libraries("${PROJECT_NAME}_RunAllTests"
        ${PROJECT_NAME}
        amolqc
        utils
        gtest gtest_main
        pthread
        )