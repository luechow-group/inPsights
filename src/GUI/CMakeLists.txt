cmake_minimum_required(VERSION 3.10)

project(GUI CXX)

set ( PROJECT_FOLDER_NAME ${PROJECT_NAME} )
set ( PROJECT_PATH ${PROJECT_NAME})
set ( OUTPUT_NAME "${PROJECT_FOLDER_NAME}" )
message ( "-- ++ CREATING PROJECT OR MAKEFILE ---> ${OUTPUT_NAME}" )


find_package(Qt5 REQUIRED COMPONENTS
        Core
        Gui
        Widgets
        3DCore
        3DExtras
        3DRender
        Charts
        )

qt5_add_resources(QT_RESOURCES myresources.qrc)

set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC ON)


set(HEADER_FILES
        include/Abstract3dObject.h
        include/Sphere.h
        include/Cylinder.h
        include/DividedCylinder.h
        include/Polyline.h
        include/Cone.h
        include/Atom3D.h
        include/Bond3D.h
        include/GuiHelper.h
        include/Electron3D.h
        include/StringMethodCoordinatesPlotter.h
        include/StringMethodValuePlotter.h
        include/AtomsVector3D.h
        include/ElectronsVector3D.h
        include/MoleculeWidget.h
        include/ParticlesVectorPath3D.h
        include/Visualization.h
        include/InPsightsWidget.h
        include/Line3D.h
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
        source/GuiHelper.cpp
        source/Electron3D.cpp
        source/AtomsVector3D.cpp
        source/ElectronsVector3D.cpp
        source/MoleculeWidget.cpp
        source/ParticlesVectorPath3D.cpp
        source/Visualization.cpp
        )

add_library(${PROJECT_NAME} STATIC
        ${SOURCE_FILES}
        ${HEADER_FILES}
        ${QT_RESOURCES}
        )
target_include_directories(${PROJECT_NAME}
        PUBLIC include
        PRIVATE source
        )

target_link_libraries(${PROJECT_NAME}
        BaseLib #includes eigen already
        Qt5::Core
        Qt5::Gui
        Qt5::Widgets
        Qt5::3DCore
        Qt5::3DExtras
        Qt5::3DRender
        Qt5::Charts
        )

add_gtest_subdirectory_with_dependency(tests "gtest;pthread")
