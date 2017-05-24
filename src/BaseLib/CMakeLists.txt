cmake_minimum_required(VERSION 3.5)

project(BaseLib CXX)


set ( PROJECT_FOLDER_NAME ${PROJECT_NAME} )

set ( PROJECT_PATH ${PROJECT_NAME})
set ( OUTPUT_NAME "${PROJECT_FOLDER_NAME}" )
message ( "-- ++ CREATING PROJECT OR MAKEFILE ---> ${OUTPUT_NAME}" )

#set(CMAKE_AUTOMOC ON)



include_directories(
        include
        source
        )

find_package(Eigen3 REQUIRED)
include_directories( ${EIGEN3_INCLUDE_DIR} )

find_package(Qt5 REQUIRED COMPONENTS
        Core
        )

set(HEADER_FILES
        include/ElementColor.h
        include/ElementInfo.h
        include/ElementTypes.h
        include/pElementInfo.h
        include/NaturalConstants.h
        #include/ChemicalSystem.h
        include/WaveFunctionParser.h
        include/Atom.h
        include/AtomCollection.h
        )

set(SOURCE_FILES
        source/pElementInfo.cpp
        source/ElementInfo.cpp
        #source/ChemicalSystem.cpp
        source/WaveFunctionParser.cpp
        source/Atom.cpp
        source/AtomCollection.cpp
        )

add_library(${PROJECT_NAME} STATIC ${SOURCE_FILES} ${HEADER_FILES} include/Atom.h)
target_link_libraries(${PROJECT_NAME}
        Qt5::Core
        #eigen
        )
target_include_directories(${PROJECT_NAME} PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        )

add_executable(${PROJECT_NAME}_main main.cpp)
target_link_libraries(${PROJECT_NAME}_main
        ${PROJECT_NAME}
        )

add_subdirectory(tests)
