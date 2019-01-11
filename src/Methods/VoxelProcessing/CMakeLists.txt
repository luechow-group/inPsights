cmake_minimum_required(VERSION 3.10)

project(VoxelProcessing CXX)

set ( PROJECT_FOLDER_NAME ${PROJECT_NAME} )
set ( PROJECT_PATH ${PROJECT_NAME})
set ( OUTPUT_NAME "${PROJECT_FOLDER_NAME}" )
message ( "-- ++ CREATING PROJECT OR MAKEFILE ---> ${OUTPUT_NAME}" )

set(HEADER_FILES
        include/Conversion.h
        include/Vertex.h
        include/Triangle.h
        include/VoxelCube.h
        include/SurfaceDataGenerator.h
        include/SurfaceData.h
        )

set(SOURCE_FILES
        source/Vertex.cpp
        source/Triangle.cpp
        source/VoxelCube.cpp
        source/SurfaceDataGenerator.cpp
        source/SurfaceData.cpp
        )

add_library(${PROJECT_NAME} STATIC ${SOURCE_FILES} ${HEADER_FILES})

target_include_directories(${PROJECT_NAME}
        PUBLIC include
        PRIVATE source
        )

target_link_libraries(${PROJECT_NAME}
        BaseLib
        dualmc
        )

add_gtest_subdirectory_with_dependency(tests "gmock;pthread")
