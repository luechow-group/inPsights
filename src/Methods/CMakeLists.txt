cmake_minimum_required(VERSION 3.6)


#============= PROJECT Methods ===================
project(Methods)

set ( PROJECT_FOLDER_NAME ${PROJECT_NAME} )

set ( PROJECT_PATH ${PROJECT_NAME})
set ( OUTPUT_NAME "${PROJECT_FOLDER_NAME}" )
message ( "-- ++ CREATING PROJECT OR MAKEFILE ---> ${OUTPUT_NAME}" )


# Find Eigen lib
find_package(Eigen3 REQUIRED)

include_directories(
        include/*
        source/*
        ${EIGEN3_INCLUDE_DIR}
        )

set(HEADER_FILES
        include/QuaternionFit.h
        include/Optimization/IOptimizer.h
        include/Optimization/LBFGSOptimizer.h
        include/Optimization/OptimizerFactory.h
        include/Optimization/NewtonRaphsonOptimizer.h
        )

set(SOURCE_FILES
        source/QuaternionFit.cpp
        source/Optimization/IOptimizer.cpp
        source/Optimization/LBFGSOptimizer.cpp
        source/Optimization/NewtonRaphsonOptimizer.cpp
        )

add_library(${PROJECT_NAME} STATIC ${SOURCE_FILES} ${HEADER_FILES})


target_include_directories(${PROJECT_NAME} PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        )

add_library(AmolqcGUI::${PROJECT_NAME} ALIAS ${PROJECT_NAME})

add_subdirectory(tests)
