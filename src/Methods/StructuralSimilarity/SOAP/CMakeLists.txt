cmake_minimum_required(VERSION 3.11)

project(SOAP CXX)

set ( PROJECT_FOLDER_NAME ${PROJECT_NAME} )
set ( PROJECT_PATH ${PROJECT_NAME})
set ( OUTPUT_NAME "${PROJECT_FOLDER_NAME}" )
message ( "-- ++ CREATING PROJECT OR MAKEFILE ---> ${OUTPUT_NAME}" )


set(HEADER_FILES
        include/RadialBasis.h
        include/AngularBasis.h
        include/NeighborhoodExpander.h
        include/Cutoff.h
        include/Environment.h
        include/SOAPSettings.h
        include/LocalSimilarity.h
        include/PowerSpectrum.h
        include/ParticleKit.h
        include/NeighborhoodExpansion.h
        include/Sinkhorn.h
        include/StructuralSimilarity.h
        include/MolecularSpectrum.h
        )

set(SOURCE_FILES
        source/RadialBasis.cpp
        source/AngularBasis.cpp
        source/Cutoff.cpp
        source/SOAPSettings.cpp
        source/ParticleKit.cpp
        source/PowerSpectrum.cpp
        source/LocalSimilarity.cpp
        source/NeighborhoodExpander.cpp
        source/NeighborhoodExpansion.cpp
        source/Sinkhorn.cpp
        source/StructuralSimilarity.cpp
        source/MolecularSpectrum.cpp
        source/Environment.cpp
        )

add_library(${PROJECT_NAME} STATIC ${SOURCE_FILES} ${HEADER_FILES})
target_include_directories(${PROJECT_NAME}
        PUBLIC include
        PRIVATE source
        )

target_link_libraries(${PROJECT_NAME}
        BaseLib
        )

add_gtest_subdirectory_with_dependency(tests "gmock;pthread")
