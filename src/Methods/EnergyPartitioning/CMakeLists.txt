cmake_minimum_required(VERSION 3.10)

project(EnergyPartitioning CXX)

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
        include/EnergyCalculator.h
        include/SimilarReferences.h
        include/OneParticleEnergies.h
        include/ClusterData.h
        #include/GlobalPermutationSorter.h
        include/GeneralStatistics.h
        include/ValueSorter.h
        include/EnergyPartitioningSettings.h
        include/VoxelCubeGeneration.h
        )

set(SOURCE_FILES
        source/RawDataReader.cpp
        source/OneParticleEnergies.cpp
        source/ClusterData.cpp
        source/GlobalClusterSorter.cpp
        source/GlobalIdentitySorter.cpp
        source/GlobalSimilaritySorter.cpp
        source/EnergyCalculator.cpp
        #source/GlobalPermutationSorter.cpp
        source/GeneralStatistics.cpp
        source/ValueSorter.cpp
        source/EnergyPartitioningSettings.cpp
        source/VoxelCubeGeneration.cpp
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
