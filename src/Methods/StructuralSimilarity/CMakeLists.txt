cmake_minimum_required(VERSION 3.5)

project(StructuralSimilarity C CXX Fortran)

set ( PROJECT_FOLDER_NAME ${PROJECT_NAME} )
set ( PROJECT_PATH ${PROJECT_NAME})
set ( OUTPUT_NAME "${PROJECT_FOLDER_NAME}" )
message ( "-- ++ CREATING PROJECT OR MAKEFILE ---> ${OUTPUT_NAME}" )


set(HEADER_FILES
        include/QuaternionFit.h
        include/Hungarian.h
        include/SimpleSorter.h
        )

set(SOURCE_FILES
        source/QuaternionFit.cpp
        source/Hungarian.cpp
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
        )

add_subdirectory(SOAP)
add_subdirectory(tests)
