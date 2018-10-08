cmake_minimum_required(VERSION 3.5)

project(LocalSpinMultiplicity C CXX Fortran)


set ( PROJECT_FOLDER_NAME ${PROJECT_NAME} )
set ( PROJECT_PATH ${PROJECT_NAME})
set ( OUTPUT_NAME "${PROJECT_FOLDER_NAME}" )
message ( "-- ++ CREATING PROJECT OR MAKEFILE ---> ${OUTPUT_NAME}" )

find_package(Eigen3 REQUIRED)
include_directories( ${EIGEN3_INCLUDE_DIR} )

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
        src/SpinDeterminer.cpp
        src/SpinProjectionQuantumNumberCounter.cpp
        )

add_library(${PROJECT_NAME} STATIC ${SOURCE_FILES} ${HEADER_FILES})
target_link_libraries(${PROJECT_NAME}
        BaseLib #includes eigen already
        )

target_include_directories(${PROJECT_NAME} PUBLIC
        PUBLIC
            include
        PRIVATE
            source
        )

add_subdirectory(tests)
