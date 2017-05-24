cmake_minimum_required(VERSION 3.5)
project(LocalSpinMultiplicity)

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR})
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(include)

find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

if(APPLE)
    add_subdirectory(../../../cmake-build-debug/googletest-src ../../../cmake-build-debug/googletest-build)
else()
        configure_file(DownloadGoogleTest.in
                googletest-download/CMakeLists.txt)
        execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
                 WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download )
        execute_process(COMMAND ${CMAKE_COMMAND} --build .
                 WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download )

# Prevent GoogleTest from overriding our compiler/linker options
# when building with Visual Studio
        set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Add googletest directly to our build. This adds
# the following targets: gtest, gtest_main, gmock
# and gmock_main

        add_subdirectory(${CMAKE_BINARY_DIR}/googletest-src
                ${CMAKE_BINARY_DIR}/googletest-build)

endif(APPLE)



#find_package(GTest REQUIRED)

set(MAIN_FILE src/main.cpp)

set(SOURCE_FILES
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
        src/ElectronPermutation.cpp include/ElectronPermutation.h include/Assignation.h src/ElectronAssigner.cpp include/ElectronAssigner.h src/HungarianElectronAssigner.cpp include/HungarianElectronAssigner.h src/TestElectronAssigner.cpp include/TestElectronAssigner.h src/Hungarian.cpp include/Hungarian.h src/SpinDeterminer.cpp include/SpinDeterminer.h src/SpinQuantumNumberCounter.cpp src/SpinQuantumNumberCounter.h)


add_library(${PROJECT_NAME}_lib STATIC ${SOURCE_FILES} ${EIGEN3_INCLUDE_DIR} ${Boost_INCLUDE_DIRS})

add_subdirectory(src/tests)

add_executable(LocalSpinMultiplicity
        ${MAIN_FILE}
        )

target_link_libraries(LocalSpinMultiplicity ${PROJECT_NAME}_lib)

#target_link_libraries(gtest gtest_main)
