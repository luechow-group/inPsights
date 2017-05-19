cmake_minimum_required(VERSION 3.5)
project(LocalSpinMultiplicity)

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR})
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(include)


set(SOURCE_FILES
        src/main.cpp
        include/pse.h
        src/Particle.cpp
        include/Particle.h
        src/Electron.cpp
        include/Electron.h
        src/Core.cpp
        include/Core.h
        src/InputOutput.cpp
        include/InputOutput.h
        src/FileXyzInput.cpp
        include/FileXyzInput.h
        src/Molecule.cpp
        include/Molecule.h
        src/pse.cpp
        src/ElectronRefArrangement.cpp
        include/ElectronRefArrangement.h
        src/ElectronPermutation.cpp include/ElectronPermutation.h include/Assignation.h src/ElectronAssigner.cpp include/ElectronAssigner.h src/HungarianElectronAssigner.cpp include/HungarianElectronAssigner.h src/TestElectronAssigner.cpp include/TestElectronAssigner.h src/Hungarian.cpp include/Hungarian.h)

add_executable(LocalSpinMultiplicity
        ${SOURCE_FILES}
        ${EIGEN3_INCLUDE_DIR}
        )
