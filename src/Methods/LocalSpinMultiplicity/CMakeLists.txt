cmake_minimum_required(VERSION 3.5)
project(LocalSpinMultiplicity)

set(CMAKE_CXX_STANDARD 11)

if ( CMAKE_COMPILER_IS_GNUCC )
        set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wno-ignored-attributes -Wno-deprecated-declarations")
endif ( CMAKE_COMPILER_IS_GNUCC )
if ( MSVC )
        set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS}")
endif ( MSVC )

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR})
find_package(Eigen3 REQUIRED)  #Eigen3 is Linear Algebra header only library
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(include)

#MESSAGE( STATUS "CMAKE_BUILD_TYPE: " ${CMAKE_BUILD_TYPE} )
STRING(TOUPPER ${CMAKE_BUILD_TYPE} CMAKE_BUILD_TYPE_UPPER)
if(CMAKE_BUILD_TYPE_UPPER MATCHES "DEBUG")
        #set(GOOGTEST_ALREADY 1) #activate this, if googletest is already available on your system. Change below folder then.
        if(GOOGTEST_ALREADY)
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

        endif(GOOGTEST_ALREADY)
endif(CMAKE_BUILD_TYPE_UPPER MATCHES "DEBUG")


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
        include/Assignment.h
        src/ElectronAssigner.cpp
        include/ElectronAssigner.h
        src/HungarianElectronAssigner.cpp
        include/HungarianElectronAssigner.h
        src/TestElectronAssigner.cpp
        include/TestElectronAssigner.h
        src/Hungarian.cpp
        include/Hungarian.h
        src/SpinDeterminer.cpp
        include/SpinDeterminer.h
        src/SpinProjectionQuantumNumberCounter.cpp
        include/SpinProjectionQuantumNumberCounter.h)


add_library(${PROJECT_NAME}_lib STATIC ${SOURCE_FILES} ${EIGEN3_INCLUDE_DIR})


if(CMAKE_BUILD_TYPE_UPPER MATCHES "DEBUG")
        add_subdirectory(src/tests)
endif(CMAKE_BUILD_TYPE_UPPER MATCHES "DEBUG")

add_executable(LocalSpinMultiplicity
        ${MAIN_FILE}
        )

target_link_libraries(LocalSpinMultiplicity ${PROJECT_NAME}_lib)

#target_link_libraries(gtest gtest_main)
