cmake_minimum_required(VERSION 3.11)

project(StructuralSimilarity CXX)

set ( PROJECT_FOLDER_NAME ${PROJECT_NAME} )
set ( PROJECT_PATH ${PROJECT_NAME})
set ( OUTPUT_NAME "${PROJECT_FOLDER_NAME}" )
message ( "-- ++ CREATING PROJECT OR MAKEFILE ---> ${OUTPUT_NAME}" )


set(HEADER_FILES
        include/QuaternionFit.h
        include/Hungarian.h
        include/SimpleSorter.h
        include/BestMatch.h
        include/GraphAnalysis.h
        include/SpinPairClassification.h
        include/Motif.h
        include/Motifs.h
        include/BestMatchSimilarity.h
        include/BestMatchDistance.h
        include/NearestElectrons.h
        )

set(SOURCE_FILES
        source/QuaternionFit.cpp
        source/SimpleSorter.cpp
        source/BestMatch.cpp
        source/GraphAnalysis.cpp
        source/SpinPairClassification.cpp
        source/Motif.cpp
        source/Motifs.cpp
        source/BestMatchSimilarity.cpp
        source/NearestElectrons.cpp
        )

add_library(${PROJECT_NAME} STATIC ${SOURCE_FILES} ${HEADER_FILES})
target_include_directories(${PROJECT_NAME}
        PUBLIC include
        PRIVATE source
        )

target_link_libraries(${PROJECT_NAME}
        BaseLib
        SOAP
        DensityBasedScan
        )

add_gtest_subdirectory_with_dependency(tests "gmock;pthread")

add_subdirectory(SOAP)
add_subdirectory(DensityBasedScan)

