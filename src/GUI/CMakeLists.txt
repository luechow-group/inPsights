cmake_minimum_required(VERSION 3.11)

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
        include/Bond3D.h
        include/GuiHelper.h
        include/MoleculeWidget.h
        include/ParticlesVectorPath3D.h
        include/Visualization.h
        include/InPsightsWidget.h
        include/Line3D.h
        include/Particle3D.h
        include/ParticlesVector3D.h
        include/IConnection.h
        include/Bonds3D.h
        include/SpinCorrelations3D.h
        include/IntegerSortedTreeWidgetItem.h
        include/SurfaceMesh.h
        include/Surface.h
        include/Arrow.h
        include/CartesianAxes.h
        include/ColorPalette.h
        include/X3domConverter.h
        )

set(SOURCE_FILES
        source/Abstract3dObject.cpp
        source/Sphere.cpp
        source/Cylinder.cpp
        source/DividedCylinder.cpp
        source/Polyline.cpp
        source/Cone.cpp
        source/Bond3D.cpp
        source/GuiHelper.cpp
        source/MoleculeWidget.cpp
        source/ParticlesVectorPath3D.cpp
        source/Visualization.cpp
        source/ParticlesVector3D.cpp
        source/Line3D.cpp
        source/Bonds3D.cpp
        source/SpinCorrelations3D.cpp
        source/InPsightsWidget.cpp
        source/IntegerSortedTreeWidgetItem.cpp
        source/SurfaceMesh.cpp
        source/Surface.cpp
        source/Arrow.cpp
        source/CartesianAxes.cpp
        source/X3domConverter.cpp
        source/ColorPalette.cpp
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
        BaseLib # includes eigen
        MaximaProcessing
        SurfaceGeneration # includes dualmc
        Qt5::Core
        Qt5::Gui
        Qt5::Widgets
        Qt5::3DCore
        Qt5::3DExtras
        Qt5::3DRender
        Qt5::Charts
        )

add_gtest_subdirectory_with_dependency(tests "gmock;pthread")
