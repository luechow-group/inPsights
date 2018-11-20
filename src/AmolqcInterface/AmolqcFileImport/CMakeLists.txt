cmake_minimum_required(VERSION 3.10)

project(AmolqcFileImport CXX)

set ( PROJECT_FOLDER_NAME ${PROJECT_NAME})
set ( PROJECT_PATH ${PROJECT_NAME})
set ( OUTPUT_NAME "${PROJECT_FOLDER_NAME}" )
message ( "-- ++ CREATING PROJECT OR MAKEFILE ---> ${OUTPUT_NAME}" )

set(HEADER_FILES
        include/WfFileImporter.h
        include/RefFileImporter.h
        include/OptimizationPathFileImporter.h
        include/AmolqcImporter.h
        )

set(SOURCE_FILES
        source/RefFileImporter.cpp
        source/WfFileImporter.cpp
        source/OptimizationPathFileImporter.cpp
        source/AmolqcImporter.cpp
        )


add_library(${PROJECT_NAME} STATIC ${SOURCE_FILES} ${HEADER_FILES})

target_include_directories(${PROJECT_NAME}
        PUBLIC
            include
        PRIVATE
            source
        )

target_link_libraries(${PROJECT_NAME}
        BaseLib
        )

#add_gtest_subdirectory_with_dependency(tests "gtest;pthread")
