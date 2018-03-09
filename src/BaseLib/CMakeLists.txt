cmake_minimum_required(VERSION 3.5)

project(BaseLib CXX)


set ( PROJECT_FOLDER_NAME ${PROJECT_NAME} )

set ( PROJECT_PATH ${PROJECT_NAME})
set ( OUTPUT_NAME "${PROJECT_FOLDER_NAME}" )
message ( "-- ++ CREATING PROJECT OR MAKEFILE ---> ${OUTPUT_NAME}" )

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
        include/pElementInfo.h
        include/ElementType.h
        include/ElementTypesVector.h
        include/NaturalConstants.h
        include/Atom.h
        include/AtomsVector.h
        include/AtomsVectorCollection.h
        include/Particle.h
        include/ParticlesVector.h
        include/ParticlesVectorCollection.h
        include/Electron.h
        include/ElectronsVector.h
        include/ElectronsVectorCollection.h
        include/SpinType.h
        include/SpinTypesVector.h
        include/Importer.h
        include/CollectionParser.h
        include/AbstractVector.h
        include/PositionsVector.h
        include/Position.h
        include/PositionsVectorCollection.h
        include/ToString.h
        )

set(SOURCE_FILES
        source/pElementInfo.cpp
        source/ElementInfo.cpp
        source/ElementTypesVector.cpp
        source/Atom.cpp
        source/AtomsVector.cpp
        source/AtomsVectorCollection.cpp
        source/Particle.cpp
        source/ParticlesVector.cpp
        source/ParticlesVectorCollection.cpp
        source/Electron.cpp
        source/ElectronsVector.cpp
        source/ElectronsVectorCollection.cpp
        source/SpinTypesVector.cpp
        source/Importer.cpp
        source/CollectionParser.cpp
        source/SpinType.cpp
        source/ElementType.cpp
        source/PositionsVector.cpp
        source/Position.cpp
        source/AbstractVector.cpp
        source/PositionsVectorCollection.cpp
        source/ToString.cpp
        )

file(GLOB_RECURSE RESOURCE_FILES resources/*.*)
file(COPY ${RESOURCE_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_library(${PROJECT_NAME} STATIC ${SOURCE_FILES} ${HEADER_FILES})
target_link_libraries(${PROJECT_NAME}
        nlohmann_json
        )
target_include_directories(${PROJECT_NAME} PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        )

add_subdirectory(tests)
