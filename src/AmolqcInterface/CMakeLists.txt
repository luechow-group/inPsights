cmake_minimum_required(VERSION 3.5)

project(AmolqcInterface C CXX Fortran)

set(Amolqcpp_BUILD TRUE)

set ( PROJECT_FOLDER_NAME ${PROJECT_NAME})
set ( PROJECT_PATH ${PROJECT_NAME})
set ( OUTPUT_NAME "${PROJECT_FOLDER_NAME}" )
message ( "-- ++ CREATING PROJECT OR MAKEFILE ---> ${OUTPUT_NAME}" )

set(CMAKE_CXX_STANDARD 11)

include_directories(
        include
        source
        )

set(HEADER_FILES
        include/ElectronicWaveFunction.h
        include/ElectronicWaveFunctionProblem.h
        include/WfFileImporter.h
        include/RefFileImporter.h
        include/OptimizationPathFileImporter.h
        include/AmolqcImporter.h
        )

set(SOURCE_FILES
        source/ElectronicWaveFunction.cpp
        source/ElectronicWaveFunctionProblem.cpp
        source/RefFileImporter.cpp
        source/WfFileImporter.cpp
        source/OptimizationPathFileImporter.cpp
        source/AmolqcImporter.cpp
        )

add_subdirectory(amolqc)

add_library(${PROJECT_NAME} STATIC ${SOURCE_FILES} ${HEADER_FILES})
target_link_libraries(${PROJECT_NAME}
        BaseLib
        cppoptlib
        AmolqcppInterface
        )
target_include_directories(${PROJECT_NAME} PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        )

add_executable(${PROJECT_NAME}_main main.cpp)
target_link_libraries(${PROJECT_NAME}_main
        ${PROJECT_NAME}
        )
