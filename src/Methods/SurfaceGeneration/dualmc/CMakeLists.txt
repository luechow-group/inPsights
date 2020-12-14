cmake_minimum_required(VERSION 3.9)
project (dualmc CXX)

set ( PROJECT_FOLDER_NAME ${PROJECT_NAME})
set ( PROJECT_PATH ${PROJECT_NAME})
set ( OUTPUT_NAME "${PROJECT_FOLDER_NAME}" )
message ( "-- ++ CREATING PROJECT OR MAKEFILE ---> ${OUTPUT_NAME}" )

set(HEADER_FILES
        include/DualMC.h
        include/Tables.h
        include/HelperFunctions.h
        )

set(SOURCE_FILES
        source/DualMC.cpp
        )

add_library(${PROJECT_NAME} STATIC ${SOURCE_FILES} ${HEADER_FILES})
target_include_directories(${PROJECT_NAME}
        PUBLIC include
        PRIVATE source
        )
