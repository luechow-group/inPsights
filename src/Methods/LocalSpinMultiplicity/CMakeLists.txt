cmake_minimum_required(VERSION 3.5)

project(LocalSpinMultiplicity C CXX Fortran)


set ( PROJECT_FOLDER_NAME ${PROJECT_NAME} )
set ( PROJECT_PATH ${PROJECT_NAME})
set ( OUTPUT_NAME "${PROJECT_FOLDER_NAME}" )
message ( "-- ++ CREATING PROJECT OR MAKEFILE ---> ${OUTPUT_NAME}" )

set(CMAKE_CXX_STANDARD 11)

find_package(Eigen3 REQUIRED)
include_directories( ${EIGEN3_INCLUDE_DIR} )

include_directories(
        include
        source
        resources
)

set(HEADER_FILES
        include/pse.h
        include/Particle.h
        include/Electron.h
        include/Core.h
        include/InputOutput.h
        include/FileXyzInput.h
        include/Molecule.h
        include/Assignment.h
        include/HungarianElectronAssigner.h
        include/Hungarian.h
        include/SpinDeterminer.h
        include/SpinProjectionQuantumNumberCounter.h
        )

set(SOURCE_FILES
        src/Particle.cpp
        src/Electron.cpp
        src/Core.cpp
        src/InputOutput.cpp
        src/FileXyzInput.cpp
        src/Molecule.cpp
        src/pse.cpp
        src/HungarianElectronAssigner.cpp
        src/Hungarian.cpp
        src/SpinDeterminer.cpp
        src/SpinProjectionQuantumNumberCounter.cpp
        )

file(GLOB_RECURSE RESOURCE_FILES resources/*.*)
file(COPY ${RESOURCE_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_library(${PROJECT_NAME} STATIC ${SOURCE_FILES} ${HEADER_FILES})
target_link_libraries(${PROJECT_NAME}
        BaseLib
        )

target_include_directories(${PROJECT_NAME} PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        )

#target_link_libraries(${PROJECT_NAME}_main
#        ${PROJECT_NAME}
#        )

add_subdirectory(tests)
