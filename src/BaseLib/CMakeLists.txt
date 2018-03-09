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
        include/ElementTypeCollection.h
        include/NaturalConstants.h
        include/Atom.h
        include/AtomCollection.h
        include/AtomsVectorCollection.h
        include/Particle.h
        include/ParticleCollection.h
        include/ParticleCollections.h
        include/Electron.h
        include/ElectronCollection.h
        include/ElectronsVectorCollection.h
        include/SpinType.h
        include/SpinTypeCollection.h
        include/Importer.h
        include/CollectionParser.h
        include/AbstractCollection.h
        include/PositionCollection.h
        include/Position.h
        include/PositionCollections.h
        include/ToString.h
        )

set(SOURCE_FILES
        source/pElementInfo.cpp
        source/ElementInfo.cpp
        source/ElementTypeCollection.cpp
        source/Atom.cpp
        source/AtomCollection.cpp
        source/AtomsVectorCollection.cpp
        source/Particle.cpp
        source/ParticleCollection.cpp
        source/ParticleCollections.cpp
        source/Electron.cpp
        source/ElectronCollection.cpp
        source/ElectronsVectorCollection.cpp
        source/SpinTypeCollection.cpp
        source/Importer.cpp
        source/CollectionParser.cpp
        source/SpinType.cpp
        source/ElementType.cpp
        source/PositionCollection.cpp
        source/Position.cpp
        source/AbstractCollection.cpp
        source/PositionCollections.cpp
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
