cmake_minimum_required(VERSION 3.5)

project(StructuralSimilarity CXX)

set ( PROJECT_FOLDER_NAME ${PROJECT_NAME} )
set ( PROJECT_PATH ${PROJECT_NAME})
set ( OUTPUT_NAME "${PROJECT_FOLDER_NAME}" )
message ( "-- ++ CREATING PROJECT OR MAKEFILE ---> ${OUTPUT_NAME}" )


set(HEADER_FILES
        include/QuaternionFit.h
        include/Hungarian.h
        include/SimpleSorter.h
        include/HungarianHelper.h
        )

set(SOURCE_FILES
        source/QuaternionFit.cpp
        source/SimpleSorter.cpp
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

add_subdirectory(SOAP)
add_subdirectory(DensityBasedScan)
add_subdirectory(tests)
