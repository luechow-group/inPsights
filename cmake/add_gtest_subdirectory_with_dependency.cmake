function(add_gtest_subdirectory_with_dependency directoryName libraries)
    FILE(GLOB TEST_SOURCE_FILES ${directoryName}/*)

    add_executable(${PROJECT_NAME}_${directoryName}.exe ${TEST_SOURCE_FILES})

    add_dependencies(inPsightsTests ${PROJECT_NAME}_${directoryName}.exe)

    target_link_libraries(${PROJECT_NAME}_${directoryName}.exe
            ${PROJECT_NAME}
            ${libraries}
            )

    add_test(NAME ${PROJECT_NAME}_${directoryName}
            COMMAND ${PROJECT_NAME}_${directoryName}.exe)

    set_tests_properties(${PROJECT_NAME}_${directoryName}
            PROPERTIES LABELS ${PROJECT_NAME})
endfunction(add_gtest_subdirectory_with_dependency)
