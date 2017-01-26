cmake_minimum_required(VERSION 3.6)
enable_testing()

find_package(GTest REQUIRED)
find_package(Eigen3 REQUIRED)

include_directories(SYSTEM
        ${GTEST_INCLUDE_DIRS}
        ${EIGEN3_INCLUDE_DIR}
        )

FILE(GLOB TEST_SOURCE_FILES unit_tests/*)

add_executable("${PROJECT_NAME}_RunAllTests" ${TEST_SOURCE_FILES})

target_link_libraries("${PROJECT_NAME}_RunAllTests"
        ${PROJECT_NAME}
        gtest gtest_main
        pthread
        )

