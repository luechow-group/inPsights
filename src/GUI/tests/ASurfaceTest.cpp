// Copyright (C) 2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <gmock/gmock.h>
#include <SurfaceDataGenerator.h>
#include <SurfaceData.h>
#include <MoleculeWidget.h>

#include <VoxelCube.h>
#include <Cylinder.h>
#include <Conversion.h>

#include <QApplication>
#include <QWidget>
#include <Qt3DExtras/Qt3DWindow>
#include <Qt3DExtras/QOrbitCameraController>
#include <Qt3DRender/QCamera>
#include <Qt3DCore/QEntity>
#include <QuickHull.h>

class RadialGaussian {
public:
    RadialGaussian(
            float variance = 1,
            float cX = 0,
            float cY = 0,
            float cZ = 0
    ) : cX(cX), cY(cY), cZ(cZ) {
        float constexpr TWO_PI = 6.283185307179586f;
        normalization = 1.0f / sqrt(TWO_PI * variance);
        falloff = -0.5f / variance;
    }

    float eval(float x, float y, float z) const {
        // compute squared input point distance to gauss center
        float const dx = x - cX;
        float const dy = y - cY;
        float const dz = z - cZ;
        float const dSquared = dx * dx + dy * dy + dz * dz;
        // compute gauss
        return normalization * exp(falloff * dSquared);
    }

private:
    // Coordinates of the sphere center.
    float cX;
    float cY;
    float cZ;
    // precomputed factors
    float normalization;
    float falloff;
};

void create3DWindow(const SurfaceData &surfaceData) {
    int argc = 0;
    char *argv[] = {};
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
    new Cylinder(entity, Qt::red, {{-0.5f, 0, 0},
                                   {0.5f,  0, 0}}, 0.005f);
    new Cylinder(entity, Qt::green, {{0, -0.5f, 0},
                                     {0, 0.5f,  0}}, 0.005f);
    new Cylinder(entity, Qt::blue, {{0, 0, -0.5f},
                                    {0, 0, 0.5f}}, 0.005f);
    new Surface(entity, surfaceData, Qt::gray);

    qt3DWindow->setRootEntity(rootEntity);
    widget->show();

    QApplication::exec();
}

TEST(ASurfaceTest, SphereFromDensity) {
    VoxelCube cube(32, 1.0, {0.5,0,0});

    RadialGaussian a(0.02,
            VoxelCube::offset_,
            VoxelCube::offset_,
            VoxelCube::offset_);
    int32_t p = 0;
    for (int32_t z = 0; z < cube.getDimensions()[2]; ++z) {
        float const nZ = float(z) * cube.inverseDimensions_[2];
        for (int32_t y = 0; y < cube.getDimensions()[1]; ++y) {
            float const nY = float(y) * cube.inverseDimensions_[1];
            for (int32_t x = 0; x < cube.getDimensions()[0]; ++x, ++p) {
                float const nX = float(x) * cube.inverseDimensions_[0];

                float rho = a.eval(nX, nY, nZ);
                if (rho > 1.0f) rho = 1.0f;

                cube.data_[p] = uint16_t(rho * std::numeric_limits<uint16_t>::max());
            }
        }
    }

    SurfaceDataGenerator surfaceDataGenerator(cube);
    auto surfaceData = surfaceDataGenerator.computeSurfaceData(0.4);

    create3DWindow(surfaceData);
}

TEST(ASurfaceTest, CubeFromData) {

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

    create3DWindow(surfaceData);
}

TEST(ASurfaceTest, CubeFromConvexHullOfPointCloud){

    std::vector<Eigen::Vector3f> pc({
        {-1,  1,  1},
        { 1,  1,  1},
        { 1,  1, -1},
        { 1, -1,  1},
        {-1, -1,  1},
        {-1, -1, -1},
        {-1,  1, -1},
        { 1, -1, -1},
        { 0,  0,  0},// Points inside the convex hull
        Eigen::Vector3f::Random(),
        Eigen::Vector3f::Random(),
        Eigen::Vector3f::Random(),
        Eigen::Vector3f::Random()
    });

    quickhull::QuickHull<float> qh;

    auto hull = qh.getConvexHull(pc);

    auto vertices = hull.getVertices();
    auto triangles = hull.getTriangles();

    Conversion::calculateVertexNormals(vertices, triangles);
    SurfaceData surfaceData;
    surfaceData.triangles = triangles;
    surfaceData.vertices = vertices;

    create3DWindow(surfaceData);
}