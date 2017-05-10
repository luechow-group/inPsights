cmake_minimum_required(VERSION 3.5)

project(GUI CXX)


set ( PROJECT_FOLDER_NAME ${PROJECT_NAME} )

set ( PROJECT_PATH ${PROJECT_NAME})
set ( OUTPUT_NAME "${PROJECT_FOLDER_NAME}" )
message ( "-- ++ CREATING PROJECT OR MAKEFILE ---> ${OUTPUT_NAME}" )

# Instruct CMake to run moc automatically when needed. Required for Qt to detect Q_OBJECT
set(CMAKE_AUTOMOC ON)


# Find Eigen lib
include_directories(
        include
        source
        )


find_package(Qt5 REQUIRED COMPONENTS
        Core
        Gui
        Widgets
        3DCore
        3DExtras
        3DRender
        #3DInput
        )

set(HEADER_FILES
        include/NaturalConstants.h
        include/ElementTypes.h
        include/ElementColor.h
        include/ElementInfo.h
        include/pElementInfo.h
        include/Abstract3dObject.h
        include/Sphere.h
        include/Cylinder.h
        include/DividedCylinder.h
        include/Atom.h
        include/Bond.h
        include/Helper.h
        include/Electron.h
        )

set(SOURCE_FILES
        source/ElementInfo.cpp
        source/pElementInfo.cpp
        source/Abstract3dObject.cpp
        source/Sphere.cpp
        source/Cylinder.cpp
        source/DividedCylinder.cpp
        source/Atom.cpp
        source/Bond.cpp
        source/Electron.cpp
        )

add_executable("${PROJECT_NAME}_main" ${SOURCE_FILES} ${HEADER_FILES} main.cpp)


# links the libraries to the executable
target_link_libraries("${PROJECT_NAME}_main"
        BSplines
        eigen
        Qt5::Core
        Qt5::Gui
        Qt5::Widgets
        Qt5::3DCore
        Qt5::3DExtras
        Qt5::3DRender
        #Qt5::3DInput
        #pthread #TODO needed?
        )
