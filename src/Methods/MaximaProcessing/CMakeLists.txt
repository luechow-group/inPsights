cmake_minimum_required(VERSION 3.11)

project(MaximaProcessing CXX)

set ( PROJECT_FOLDER_NAME ${PROJECT_NAME} )
set ( PROJECT_PATH ${PROJECT_NAME})
set ( OUTPUT_NAME "${PROJECT_FOLDER_NAME}" )
message ( "-- ++ CREATING PROJECT OR MAKEFILE ---> ${OUTPUT_NAME}" )

set(HEADER_FILES
        include/RawDataReader.h
        include/Sample.h
        include/Maximum.h
        include/ClusterData.h
        include/IdentityClusterer.h
        include/PreClusterer.h
        include/DensityBasedClusterer.h
        include/SOAPClusterer.h
        include/MaximaProcessor.h
        include/GeneralStatistics.h
        include/MaximaProcessingSettings.h
        include/VoxelCubeGeneration.h
        include/EnergyPartitioning.h
        include/EnergyStatistics.h
        include/Cluster.h
        include/IProcess.h
        include/ReferencePositionsClusterer.h
        include/VoxelCubeOverlapCalculation.h
        include/SpinCorrelationValueHistogram.h
        include/ClusterNumberAnalyzer.h
        include/ClusterAnalysis.h
        include/TotalWeightDifferenceAnalyzer.h
        include/SelectionEnergyCalculator.h
        include/EnergyResultsBundle.h
        include/ClusteringMetric.h
        )

set(SOURCE_FILES
        source/RawDataReader.cpp
        source/Sample.cpp
        source/Maximum.cpp
        source/ClusterData.cpp
        source/IdentityClusterer.cpp
        source/PreClusterer.cpp
        source/DensityBasedClusterer.cpp
        source/SOAPClusterer.cpp
        source/MaximaProcessor.cpp
        source/GeneralStatistics.cpp
        source/MaximaProcessingSettings.cpp
        source/VoxelCubeGeneration.cpp
        source/EnergyPartitioning.cpp
        source/EnergyStatistics.cpp
        source/Cluster.cpp
        source/IProcess.cpp
        source/ReferencePositionsClusterer.cpp
        source/VoxelCubeOverlapCalculation.cpp
        source/SpinCorrelationValueHistogram.cpp
        source/ClusterNumberAnalyzer.cpp
        source/ClusterAnalysis.cpp
        source/TotalWeightDifferenceAnalyzer.cpp
        source/SelectionEnergyCalculator.cpp
        source/ClusteringMetric.cpp
        )

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
