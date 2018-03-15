cmake_minimum_required(VERSION 3.5)

project(Executables C CXX Fortran)


set ( PROJECT_FOLDER_NAME ${PROJECT_NAME} )

set ( PROJECT_PATH ${PROJECT_NAME})
set ( OUTPUT_NAME "${PROJECT_FOLDER_NAME}" )
message ( "-- ++ CREATING PROJECT OR MAKEFILE ---> ${OUTPUT_NAME}" )

# Instruct CMake to run moc automatically when needed. Required for Qt to detect Q_OBJECT
set(CMAKE_AUTOMOC ON)

# Find Eigen lib
find_package(Eigen3 REQUIRED)
include_directories( ${EIGEN3_INCLUDE_DIR} )

include_directories(
        resources
        )

file(GLOB_RECURSE RESOURCE_FILES resources/*.*)
file(COPY ${RESOURCE_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_executable("${PROJECT_NAME}_mainStringMethod" ${SOURCE_FILES} ${HEADER_FILES} mainStringMethodGUI.cpp)
add_executable("${PROJECT_NAME}_mainPathOptimization" ${SOURCE_FILES} ${HEADER_FILES} mainPathOptimization.cpp)
add_executable("${PROJECT_NAME}_mainHessianTesting" ${SOURCE_FILES} ${HEADER_FILES} mainHessianTesting.cpp)
add_executable("${PROJECT_NAME}_mainLocalPotentialEnergy" ${SOURCE_FILES} ${HEADER_FILES} mainLocalPotentialEnergy.cpp)

add_library(${PROJECT_NAME} INTERFACE
        )

target_link_libraries(${PROJECT_NAME} INTERFACE
        GUI
        BaseLib
        AmolqcInterface
        LocalPotentialEnergy
        TransitionStateSearch
        )

target_include_directories(${PROJECT_NAME} INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/resources>
        )

# links the libraries to the executable
target_link_libraries("${PROJECT_NAME}_mainStringMethod"
        ${PROJECT_NAME}
        )

target_link_libraries("${PROJECT_NAME}_mainPathOptimization"
        ${PROJECT_NAME}
        )

target_link_libraries("${PROJECT_NAME}_mainHessianTesting"
        ${PROJECT_NAME}
        )

target_link_libraries("${PROJECT_NAME}_mainLocalPotentialEnergy"
        ${PROJECT_NAME}
        )



