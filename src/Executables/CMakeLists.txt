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

add_executable("StringMethodGUI.exe" ${SOURCE_FILES} ${HEADER_FILES} StringMethodGUI.cpp)
add_executable("PathOptimization.exe" ${SOURCE_FILES} ${HEADER_FILES} PathOptimization.cpp)
add_executable("HessianTesting.exe" ${SOURCE_FILES} ${HEADER_FILES} HessianTesting.cpp)
add_executable("LocalPotentialEnergy.exe" ${SOURCE_FILES} ${HEADER_FILES} LocalPotentialEnergy.cpp)

add_library(Resources INTERFACE
        )

target_include_directories(Resources INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/resources>
        )

# links the libraries to the executable
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



