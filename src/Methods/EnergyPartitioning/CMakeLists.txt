cmake_minimum_required(VERSION 3.5)

project(EnergyPartitioning CXX)
set(CMAKE_CXX_STANDARD 17)

set ( PROJECT_FOLDER_NAME ${PROJECT_NAME} )
set ( PROJECT_PATH ${PROJECT_NAME})
set ( OUTPUT_NAME "${PROJECT_FOLDER_NAME}" )
message ( "-- ++ CREATING PROJECT OR MAKEFILE ---> ${OUTPUT_NAME}" )


set(HEADER_FILES
        include/RawDataReader.h
        include/SampleData.h
        include/ReferenceData.h
        )

set(SOURCE_FILES
        source/RawDataReader.cpp
        )

add_library(${PROJECT_NAME} STATIC ${SOURCE_FILES} ${HEADER_FILES})
target_include_directories(${PROJECT_NAME}
        PUBLIC include
        PRIVATE source
        )

target_link_libraries(${PROJECT_NAME}
        BaseLib
        StructuralSimilarity
        )

add_subdirectory(tests)
