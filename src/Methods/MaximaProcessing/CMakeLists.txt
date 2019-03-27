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
        include/GlobalIdentitySorter.h
        include/GlobalSimilaritySorter.h
        include/GlobalClusterSorter.h
        include/MaximaProcessor.h
        include/SimilarReferences.h
        include/ClusterData.h
        #include/GlobalPermutationSorter.h
        include/GeneralStatistics.h
        include/ValueSorter.h
        include/MaximaProcessingSettings.h
        include/VoxelCubeGeneration.h
        include/EnergyPartitioning.h
        include/EnergyStatistics.h
        )

set(SOURCE_FILES
        source/RawDataReader.cpp
        source/ClusterData.cpp
        source/GlobalClusterSorter.cpp
        source/GlobalIdentitySorter.cpp
        source/GlobalSimilaritySorter.cpp
        source/MaximaProcessor.cpp
        #source/GlobalPermutationSorter.cpp
        source/GeneralStatistics.cpp
        source/ValueSorter.cpp
        source/MaximaProcessingSettings.cpp
        source/VoxelCubeGeneration.cpp
        source/EnergyPartitioning.cpp
        source/EnergyStatistics.cpp
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
