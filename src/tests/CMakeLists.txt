cmake_minimum_required(VERSION 3.5)
enable_testing()


#FILE(GLOB TEST_SOURCE_FILES
#        *
        #integration_tests/*
        #functional_tests/*
#        )
set(TEST_SOURCE_FILES testmain.cpp ATestElectronAssignerTest.cpp AFindElementTest.cpp AHungarianTest.cpp AFileXyzInputTest.cpp)

add_executable("${PROJECT_NAME}_RunAllTests" ${TEST_SOURCE_FILES})

target_link_libraries("${PROJECT_NAME}_RunAllTests"
        gtest gtest_main gmock
        pthread
        ${PROJECT_NAME}_lib
        )

