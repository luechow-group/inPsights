cmake_minimum_required(VERSION 3.5)

project(AmolqcInterface C CXX Fortran)

set ( PROJECT_FOLDER_NAME ${PROJECT_NAME})
set ( PROJECT_PATH ${PROJECT_NAME})
set ( OUTPUT_NAME "${PROJECT_FOLDER_NAME}" )
message ( "-- ++ CREATING PROJECT OR MAKEFILE ---> ${OUTPUT_NAME}" )

set(HEADER_FILES
        include/ElectronicWaveFunction.h
        include/ElectronicWaveFunctionProblem.h
        include/AmolqcFileImport/WfFileImporter.h
        include/AmolqcFileImport/RefFileImporter.h
        include/AmolqcFileImport/OptimizationPathFileImporter.h
        include/AmolqcFileImport/AmolqcImporter.h
        )

set(SOURCE_FILES
        source/ElectronicWaveFunction.cpp
        source/ElectronicWaveFunctionProblem.cpp
        source/AmolqcFileImport/RefFileImporter.cpp
        source/AmolqcFileImport/WfFileImporter.cpp
        source/AmolqcFileImport/OptimizationPathFileImporter.cpp
        source/AmolqcFileImport/AmolqcImporter.cpp
        )

add_subdirectory(amolqc)

add_library(${PROJECT_NAME} STATIC ${SOURCE_FILES} ${HEADER_FILES})

target_include_directories(${PROJECT_NAME}
        PUBLIC include
        PRIVATE source
        )

target_link_libraries(${PROJECT_NAME}
        BaseLib
        cppoptlib
        AmolqcppInterface
        )

add_subdirectory(tests)
