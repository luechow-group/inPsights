cmake_minimum_required(VERSION 3.5)

project(StructuralSimilarity C CXX Fortran)



set ( PROJECT_FOLDER_NAME ${PROJECT_NAME} )
set ( PROJECT_PATH ${PROJECT_NAME})
set ( OUTPUT_NAME "${PROJECT_FOLDER_NAME}" )
message ( "-- ++ CREATING PROJECT OR MAKEFILE ---> ${OUTPUT_NAME}" )

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
        include/SimpleSorter.h
        )

set(SOURCE_FILES
        source/RadialBasis.cpp
        source/AtomicNeighborhoodDensityFunction.cpp
        source/QuaternionFit.cpp
        source/Hungarian.cpp
        )

add_library(${PROJECT_NAME} STATIC ${SOURCE_FILES} ${HEADER_FILES})
target_include_directories(${PROJECT_NAME} PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        )

target_link_libraries(${PROJECT_NAME}
        BaseLib
        SOAP
        SphericalIntegration
        )

add_subdirectory(SOAP)
add_subdirectory(tests)


