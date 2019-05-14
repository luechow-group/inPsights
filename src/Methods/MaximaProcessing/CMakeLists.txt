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
        include/BestMatchDistanceIdentityClusterer.h
        include/BestMatchDistanceSimilarityClusterer.h
        include/BestMatchDistanceDensityBasedClusterer.h
        include/BestMatchSOAPSimilarityClusterer.h
        include/MaximaProcessor.h
        include/GeneralStatistics.h
        include/ValueSorter.h
        include/MaximaProcessingSettings.h
        include/VoxelCubeGeneration.h
        include/EnergyPartitioning.h
        include/EnergyStatistics.h
        include/Group.h
        include/IClusterer.h
        )

set(SOURCE_FILES
        source/RawDataReader.cpp
        source/Sample.cpp
        source/Reference.cpp
        source/ClusterData.cpp
        source/BestMatchDistanceIdentityClusterer.cpp
        source/BestMatchDistanceSimilarityClusterer.cpp
        source/BestMatchDistanceDensityBasedClusterer.cpp
        source/BestMatchSOAPSimilarityClusterer.cpp
        source/MaximaProcessor.cpp
        source/GeneralStatistics.cpp
        source/MaximaProcessingSettings.cpp
        source/VoxelCubeGeneration.cpp
        source/EnergyPartitioning.cpp
        source/EnergyStatistics.cpp
        source/Group.cpp
        source/IClusterer.cpp
        )

add_library(${PROJECT_NAME} STATIC ${SOURCE_FILES} ${HEADER_FILES})
target_include_directories(${PROJECT_NAME}
        PUBLIC include
        PRIVATE source
        )

target_link_libraries(${PROJECT_NAME}
        BaseLib
        StructuralSimilarity
        VoxelProcessing
        )

add_gtest_subdirectory_with_dependency(tests "gmock;pthread")
