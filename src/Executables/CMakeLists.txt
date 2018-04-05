cmake_minimum_required(VERSION 3.5)

project(Executables C CXX Fortran)

set ( PROJECT_FOLDER_NAME ${PROJECT_NAME} )

set ( PROJECT_PATH ${PROJECT_NAME})
set ( OUTPUT_NAME "${PROJECT_FOLDER_NAME}" )
message ( "-- ++ CREATING PROJECT OR MAKEFILE ---> ${OUTPUT_NAME}" )

include_directories(
        resources
        )

file(GLOB_RECURSE RESOURCE_FILES resources/*.*)
file(COPY ${RESOURCE_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_executable("StringMethodGUI.exe" StringMethodGUI.cpp)
add_executable("PathOptimization.exe" PathOptimization.cpp)
add_executable("HessianTesting.exe" HessianTesting.cpp)
add_executable("LocalPotentialEnergy.exe" LocalPotentialEnergy.cpp)

add_library(Resources INTERFACE
        )

target_include_directories(Resources INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/resources>
        )

target_link_libraries("StringMethodGUI.exe"
        Resources
        GUI
        BaseLib
        BSplines
        TransitionStateSearch
        )

target_link_libraries("PathOptimization.exe"
        Resources
        GUI
        BaseLib
        BSplines
        TransitionStateSearch
        )

target_link_libraries("HessianTesting.exe"
        Resources
        GUI
        BaseLib
        BSplines
        TransitionStateSearch
        )

target_link_libraries("LocalPotentialEnergy.exe"
        Resources
        GUI
        LocalPotentialEnergy
        AmolqcInterface
        )
