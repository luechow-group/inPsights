cmake_minimum_required(VERSION 3.5)


project(SphericalHarmonics CXX)

set ( PROJECT_FOLDER_NAME ${PROJECT_NAME})

set ( PROJECT_PATH ${PROJECT_NAME})
set ( OUTPUT_NAME "${PROJECT_FOLDER_NAME}" )
message ( "-- ++ CREATING PROJECT OR MAKEFILE ---> ${OUTPUT_NAME}" )

set(CMAKE_CXX_STANDARD 14)
add_definitions(
        -O3
        -Wall
        -pedantic
        -pthread
        )

set(HEADER_FILES
        include/BoostSphericalHarmonics.h
        )
set(SOURCE_FILES
        source/BoostSphericalHarmonics.cpp
        )

add_library(${PROJECT_NAME} STATIC ${HEADER_FILES} ${SOURCE_FILES})

#TODO put put library somewhere else cause it does not belong in Externals
target_link_libraries(${PROJECT_NAME} eigen)

target_include_directories(${PROJECT_NAME} PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        )

add_subdirectory(tests)
