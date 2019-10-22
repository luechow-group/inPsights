function(add_gtest_subdirectory_with_dependency directoryName libraries)
    FILE(GLOB TEST_SOURCE_FILES ${directoryName}/*)

    set(TEST_EXECUTABLE_NAME "${PROJECT_NAME}_${directoryName}")

    add_executable(${TEST_EXECUTABLE_NAME} ${TEST_SOURCE_FILES})

    add_dependencies(inPsightsTests ${TEST_EXECUTABLE_NAME})

    target_link_libraries(${TEST_EXECUTABLE_NAME}
            ${PROJECT_NAME}
            ${libraries}
            )

    add_test(NAME ${TEST_EXECUTABLE_NAME}
            COMMAND ${TEST_EXECUTABLE_NAME})

    set_tests_properties(${TEST_EXECUTABLE_NAME}
            PROPERTIES LABELS ${PROJECT_NAME})
endfunction(add_gtest_subdirectory_with_dependency)
