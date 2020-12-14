cmake_minimum_required(VERSION 3.11)

project(SurfaceGeneration CXX)

add_subdirectory(dualmc)

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

        include/VertexDataSource.h
        include/Types.h
        include/Ray.h
        include/QuickHull.h
        include/Pool.h
        include/Plane.h
        include/Mesh.h
        include/MathUtils.h
        include/HalfEdgeMesh.h
        include/ConvexHull.h
        )

set(SOURCE_FILES
        source/Conversion.cpp
        source/Vertex.cpp
        source/Triangle.cpp
        source/VoxelCube.cpp
        source/SurfaceDataGenerator.cpp
        source/SurfaceData.cpp

        source/QuickHull.cpp
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
