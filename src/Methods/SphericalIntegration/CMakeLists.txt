cmake_minimum_required(VERSION 3.5)

project(SphericalIntegration CXX)


set ( PROJECT_FOLDER_NAME ${PROJECT_NAME} )
set ( PROJECT_PATH ${PROJECT_NAME})
set ( OUTPUT_NAME "${PROJECT_FOLDER_NAME}" )
message ( "-- ++ CREATING PROJECT OR MAKEFILE ---> ${OUTPUT_NAME}" )

include_directories(
        include
        source
        )

set(HEADER_FILES
        include/SphericalIntegrator.h
        include/SphericalIntegratorSettings.h
        )

set(SOURCE_FILES
        source/SphericalIntegrator.cpp
        source/SphericalIntegratorSettings.cpp
        )

add_library(${PROJECT_NAME} STATIC ${SOURCE_FILES} ${HEADER_FILES})
target_include_directories(${PROJECT_NAME} PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        )

target_link_libraries(${PROJECT_NAME}
        BaseLib #includes eigen already
        GaussKronrodCartesianIntegration
        LebedevSphericalIntegration
        SphericalHarmonics
        )

add_subdirectory(tests)
