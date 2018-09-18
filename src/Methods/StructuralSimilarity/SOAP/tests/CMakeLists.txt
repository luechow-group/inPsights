cmake_minimum_required(VERSION 3.5)
enable_testing()

FILE(GLOB TEST_SOURCE_FILES
        *
        unit/*
        integration/*
        )

include_directories(
        ${gtest_SOURCE_DIR}/include
        ${gtest_SOURCE_DIR}
        )

add_executable("${PROJECT_NAME}Tests.exe" ${TEST_SOURCE_FILES} main.cpp)

target_link_libraries("${PROJECT_NAME}Tests.exe"
        ${PROJECT_NAME}
        gmock
        gtest
        pthread
        )