function(add_gtests)
    FILE(GLOB TEST_SOURCE_FILES *)

    add_executable(${PROJECT_NAME}Tests.exe ${TEST_SOURCE_FILES})

    add_dependencies(AmolqcppTests ${PROJECT_NAME}Tests.exe)

    target_link_libraries(${PROJECT_NAME}Tests.exe
            ${PROJECT_NAME}
            gtest
            gmock
            pthread
            )

    add_test(NAME ${PROJECT_NAME}Tests
            COMMAND ${PROJECT_NAME}Tests.exe)

    set_tests_properties(${PROJECT_NAME}Tests
            PROPERTIES LABELS ${PROJECT_NAME})
endfunction(add_gtests)
