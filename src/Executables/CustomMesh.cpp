#include <QApplication>
#include <QWidget>
#include <Qt3DExtras/Qt3DWindow>
#include <Qt3DExtras/QOrbitCameraController>
#include <Qt3DRender/QCamera>
#include <Qt3DCore/QEntity>
#include <Qt3DCore/QTransform>
#include <Qt3DRender/QGeometryRenderer>
#include <Qt3DRender/QAttribute>
#include <Qt3DRender/QBuffer>
#include <Qt3DExtras/QPhongAlphaMaterial>

#include <iostream>
#include <Surface.h>

#include <DualMC.h>
#include <Conversion.h>
#include <SurfaceData.h>
#include <Vertex.h>
#include <Triangle.h>

int main(int argc, char* argv[])
{
    std::vector<Vertex> vertices({
        {{1.0f, 1.0f, 0.0f}},
        {{0.0f, 1.0f, 0.0f}},
        {{0.0f, 0.0f, 0.0f}},
        {{1.0f, 0.0f, 0.0f}},
        {{1.0f, 0.0f, 1.0f}},
        {{0.0f, 0.0f, 1.0f}},
        {{0.0f, 1.0f, 1.0f}},
        {{1.0f, 1.0f, 1.0f}}
    });

    std::vector<dualmc::Quad> quads({
        {0, 3, 2, 1},
        {4, 5, 2, 3},
        {6, 1, 2, 5},
        {7, 4, 3, 0},
        {7, 0, 1, 6},
        {7, 6, 5, 4}
    });

    auto triangles = Conversion::quadsToTriangles(quads);
    Conversion::calculateVertexNormals(vertices, triangles);
    SurfaceData surfaceData;
    surfaceData.triangles = triangles;
    surfaceData.vertices = vertices;

    QApplication app(argc, argv);

    // Root entity
    auto *rootEntity = new Qt3DCore::QEntity();

    // Window container
    auto qt3DWindow = new Qt3DExtras::Qt3DWindow();
    qt3DWindow->setRootEntity(rootEntity);
    auto widget = QWidget::createWindowContainer(qt3DWindow);

    // Camera
    auto *camController = new Qt3DExtras::QOrbitCameraController(rootEntity);

    qt3DWindow->setRootEntity(rootEntity);
    qt3DWindow->camera()->lens()->setPerspectiveProjection(45.0f, 16.0f / 9.0f, 0.1f, 100.0f);
    qt3DWindow->camera()->setPosition(QVector3D(2.5, -8, 0.0));
    qt3DWindow->camera()->setViewCenter(QVector3D(0, 0, 0));

    // For camera controls
    camController->setLinearSpeed(50.f);
    camController->setLookSpeed(180.f);
    camController->setCamera(qt3DWindow->camera());

    // Isosurface
    auto *entity = new Qt3DCore::QEntity(rootEntity);
    new Surface(entity, surfaceData, Qt::red);

    qt3DWindow->setRootEntity(rootEntity);
    widget->show();

    return QApplication::exec();
}
