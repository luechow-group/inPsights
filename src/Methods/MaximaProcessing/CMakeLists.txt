cmake_minimum_required(VERSION 3.11)

project(MaximaProcessing CXX)

set ( PROJECT_FOLDER_NAME ${PROJECT_NAME} )
set ( PROJECT_PATH ${PROJECT_NAME})
set ( OUTPUT_NAME "${PROJECT_FOLDER_NAME}" )
message ( "-- ++ CREATING PROJECT OR MAKEFILE ---> ${OUTPUT_NAME}" )

set(HEADER_FILES
        include/RawDataReader.h
        include/Sample.h
        include/Reference.h
        include/ClusterData.h
        include/IdentityClusterer.h
        include/DistanceClusterer.h
        include/DensityBasedClusterer.h
        include/SOAPClusterer.h
        include/MaximaProcessor.h
        include/GeneralStatistics.h
        include/ValueSorter.h
        include/MaximaProcessingSettings.h
        include/VoxelCubeGeneration.h
        include/EnergyPartitioning.h
        include/EnergyStatistics.h
        include/Group.h
        include/IClusterer.h
        include/ReferencePositionsClusterer.h
        include/VoxelCubeOverlapCalculation.h
        include/SpinCorrelationValueHistogram.h include/GraphClusterer.h)

set(SOURCE_FILES
        source/RawDataReader.cpp
        source/Sample.cpp
        source/Reference.cpp
        source/ClusterData.cpp
        source/IdentityClusterer.cpp
        source/DistanceClusterer.cpp
        source/DensityBasedClusterer.cpp
        source/SOAPClusterer.cpp
        source/MaximaProcessor.cpp
        source/GeneralStatistics.cpp
        source/MaximaProcessingSettings.cpp
        source/VoxelCubeGeneration.cpp
        source/EnergyPartitioning.cpp
        source/EnergyStatistics.cpp
        source/Group.cpp
        source/IClusterer.cpp
        source/ReferencePositionsClusterer.cpp
        source/VoxelCubeOverlapCalculation.cpp
        source/SpinCorrelationValueHistogram.cpp
        source/GraphClusterer.cpp)

add_library(${PROJECT_NAME} STATIC ${SOURCE_FILES} ${HEADER_FILES})
target_include_directories(${PROJECT_NAME}
        PUBLIC include
        PRIVATE source
        )

target_link_libraries(${PROJECT_NAME}
        BaseLib
        StructuralSimilarity
        SurfaceGeneration
        )

add_gtest_subdirectory_with_dependency(tests "gmock;pthread")
