cmake_minimum_required(VERSION 3.6)


#============= PROJECT alglib ===================
project(alglib)

set ( PROJECT_FOLDER_NAME ${PROJECT_NAME})

set ( PROJECT_PATH ${PROJECT_NAME})
set ( OUTPUT_NAME "${PROJECT_FOLDER_NAME}" )
message ( "-- ++ CREATING PROJECT OR MAKEFILE ---> ${OUTPUT_NAME}" )

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-g -Wall")


include_directories(
        include
        source
        )

file(GLOB HEADER_FILES "include/*.h")

file(GLOB SOURCE_FILES "source/*.cpp")

add_library(${PROJECT_NAME} STATIC ${SOURCE_FILES} ${HEADER_FILES})

target_include_directories(${PROJECT_NAME} PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        )

add_library(Externals::${PROJECT_NAME} ALIAS ${PROJECT_NAME})

