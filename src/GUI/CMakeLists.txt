cmake_minimum_required(VERSION 3.5)

project(GUI C CXX Fortran)


set ( PROJECT_FOLDER_NAME ${PROJECT_NAME} )

set ( PROJECT_PATH ${PROJECT_NAME})
set ( OUTPUT_NAME "${PROJECT_FOLDER_NAME}" )
message ( "-- ++ CREATING PROJECT OR MAKEFILE ---> ${OUTPUT_NAME}" )

# Instruct CMake to run moc automatically when needed. Required for Qt to detect Q_OBJECT
set(CMAKE_AUTOMOC ON)


# Find Eigen lib
find_package(Eigen3 REQUIRED)
include_directories( ${EIGEN3_INCLUDE_DIR} )

include_directories(
        include
        source
        resources
        )


find_package(Qt5 REQUIRED COMPONENTS
        Core
        Gui
        Widgets
        3DCore
        3DExtras
        3DRender
        Charts
        #3DInput
        )

set(HEADER_FILES
        include/Abstract3dObject.h
        include/Sphere.h
        include/Cylinder.h
        include/DividedCylinder.h
        include/Polyline.h
        include/Cone.h
        include/Atom3D.h
        include/Bond3D.h
        include/Helper.h
        include/Electron3D.h
        include/StringMethodCoordinatesPlotter.h
        include/StringMethodValuePlotter.h
        include/AtomsVector3D.h
        include/ElectronsVector3D.h
        include/MoleculeWidget.h
        include/ParticlesVectorPath3D.h
        )

set(SOURCE_FILES
        source/Abstract3dObject.cpp
        source/Sphere.cpp
        source/Cylinder.cpp
        source/DividedCylinder.cpp
        source/Polyline.cpp
        source/Cone.cpp
        source/Atom3D.cpp
        source/Bond3D.cpp
        source/Electron3D.cpp
        source/AtomsVector3D.cpp
        source/ElectronsVector3D.cpp
        source/MoleculeWidget.cpp
        source/ParticlesVectorPath3D.cpp
        )

file(GLOB_RECURSE RESOURCE_FILES resources/*.*)
file(COPY ${RESOURCE_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_executable("${PROJECT_NAME}_mainStringMethod" ${SOURCE_FILES} ${HEADER_FILES} mainStringMethodGUI.cpp)
add_executable("${PROJECT_NAME}_mainPathOptimization" ${SOURCE_FILES} ${HEADER_FILES} mainPathOptimization.cpp)
add_executable("${PROJECT_NAME}_mainHessianTesting" ${SOURCE_FILES} ${HEADER_FILES} mainHessianTesting.cpp)
add_executable("${PROJECT_NAME}_mainLocalPotentialEnergy" ${SOURCE_FILES} ${HEADER_FILES} mainLocalPotentialEnergy.cpp)

add_library(GUILib INTERFACE
        )

target_link_libraries(GUILib INTERFACE
        BaseLib
        Qt5::Core
        Qt5::Gui
        Qt5::Widgets
        Qt5::3DCore
        Qt5::3DExtras
        Qt5::3DRender
        Qt5::Charts
        )

# links the libraries to the executable
target_link_libraries("${PROJECT_NAME}_mainStringMethod"
        GUILib
        BSplines
        TransitionStateSearch
        )

target_link_libraries("${PROJECT_NAME}_mainPathOptimization"
        GUILib
        BSplines
        TransitionStateSearch
        )

target_link_libraries("${PROJECT_NAME}_mainHessianTesting"
        GUILib
        BSplines
        TransitionStateSearch
        )

target_link_libraries("${PROJECT_NAME}_mainLocalPotentialEnergy"
        GUILib
        LocalPotentialEnergy
        )
