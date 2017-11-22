cmake_minimum_required(VERSION 3.5)
enable_testing()

FILE(GLOB TEST_SOURCE_FILES *)

include_directories(
        resources
)

file(GLOB_RECURSE RESOURCE_FILES resources/*.*)
file(COPY ${RESOURCE_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_executable("${PROJECT_NAME}_RunAllTests" ${TEST_SOURCE_FILES})

target_link_libraries("${PROJECT_NAME}_RunAllTests"
        ${PROJECT_NAME}
        gmock gmock_main
        gtest gtest_main
        pthread
        )
