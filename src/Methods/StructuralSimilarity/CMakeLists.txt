cmake_minimum_required(VERSION 3.5)

project(StructuralSimilarity C CXX Fortran)

set ( PROJECT_FOLDER_NAME ${PROJECT_NAME} )
set ( PROJECT_PATH ${PROJECT_NAME})
set ( OUTPUT_NAME "${PROJECT_FOLDER_NAME}" )
message ( "-- ++ CREATING PROJECT OR MAKEFILE ---> ${OUTPUT_NAME}" )

find_package(Eigen3 REQUIRED)
include_directories( ${EIGEN3_INCLUDE_DIR} )

if(DEFINED ENV{BOOST_ROOT})
        set(BOOST_ROOT "$ENV{BOOST_ROOT}")
        message("-- BOOST_ROOT is set: ${BOOST_ROOT}")
else(DEFINED ENV{BOOST_ROOT})
        message("-- Note: BOOST_ROOT not set. Let's see whether we can find BOOST anyway.")
endif(DEFINED ENV{BOOST_ROOT})
find_package(Boost)
include_directories(${Boost_INCLUDE_DIRS})

include_directories(
        include
        source
        resources
        )

set(HEADER_FILES
        include/RadialBasis.h
        include/AtomicNeighborhoodDensityFunction.h
        include/SphericalHarmonicsRadialBasisExpander.h
        include/QuaternionFit.h
        include/Hungarian.h
        include/RadialGaussianBasis.h
        include/Gaussian.h
        include/AngularBasis.h
        include/NeighborhoodExpander.h
        include/CutoffFunction.h
        include/Environment.h
        include/ExpansionSettings.h
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
        source/AtomicNeighborhoodDensityFunction.cpp
        source/QuaternionFit.cpp
        source/Hungarian.cpp
        source/Gaussian.cpp
        source/RadialGaussianBasis.cpp
        source/AngularBasis.cpp
        source/CutoffFunction.cpp
        source/ExpansionSettings.cpp
        source/ParticleKit.cpp
        source/PowerSpectrum.cpp
        source/LocalSimilarity.cpp
        source/NeighborhoodExpander.cpp
        source/NeighborhoodExpansion.cpp
        source/Sinkhorn.cpp
        source/StructuralSimilarity.cpp
        source/MolecularSpectrum.cpp
        )

add_library(${PROJECT_NAME} STATIC ${SOURCE_FILES} ${HEADER_FILES})
target_include_directories(${PROJECT_NAME} PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        )

target_link_libraries(${PROJECT_NAME}
        BaseLib
        SphericalIntegration
        SphericalHarmonics
        )

add_subdirectory(tests)


