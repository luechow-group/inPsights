cmake_minimum_required(VERSION 3.5)
enable_testing()

FILE(GLOB TEST_SOURCE_FILES
        unit_tests/*
        functional_tests/*
        )

add_executable("${PROJECT_NAME}_RunAllTests" ${TEST_SOURCE_FILES})

target_link_libraries("${PROJECT_NAME}_RunAllTests"
        ${PROJECT_NAME}
        gtest gtest_main
        pthread
        )

