//
// Created by Michael Heuer on 2019-01-09.
//

#include <gmock/gmock.h>
#include <SurfaceDataGenerator.h>
#include <SurfaceData.h>
#include <MoleculeWidget.h>

#include <VoxelCube.h>
#include <Cylinder.h>

#include <QApplication>
#include <QWidget>
#include <Qt3DExtras/Qt3DWindow>
#include <Qt3DExtras/QOrbitCameraController>
#include <Qt3DRender/QCamera>
#include <Qt3DCore/QEntity>

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


TEST(ASurfaceTest, SphereFromDensity) {
    VoxelCube cube(32, 1.0, {0.5,0,0});

    RadialGaussian a(0.02,
            VoxelCube::offset,
            VoxelCube::offset,
            VoxelCube::offset);
    int32_t p = 0;
    for (int32_t z = 0; z < cube.dimension; ++z) {
        float const nZ = float(z) * cube.inverseDimension;
        for (int32_t y = 0; y < cube.dimension; ++y) {
            float const nY = float(y) * cube.inverseDimension;
            for (int32_t x = 0; x < cube.dimension; ++x, ++p) {
                float const nX = float(x) * cube.inverseDimension;

                float rho = a.eval(nX, nY, nZ);
                if (rho > 1.0f) rho = 1.0f;

                cube.data[p] = uint16_t(rho * std::numeric_limits<uint16_t>::max());
            }
        }
    }

    SurfaceDataGenerator surfaceDataGenerator(cube);
    auto surfaceData = surfaceDataGenerator.computeSurfaceData(0.4);

    int argc = 0;
    char* argv[] = {};
    QApplication app(argc,argv);

    // Root entity
    auto *rootEntity = new Qt3DCore::QEntity();

    // Window container
    auto qt3DWindow = new Qt3DExtras::Qt3DWindow();
    qt3DWindow->setRootEntity(rootEntity);
    auto widget = QWidget::createWindowContainer(qt3DWindow);

    // Camera
    auto *camController = new Qt3DExtras::QOrbitCameraController(rootEntity);
    qt3DWindow->camera()->lens()->setPerspectiveProjection(45.0f, 16.0f / 9.0f, 0.1f, 100.0f);
    qt3DWindow->camera()->setPosition(QVector3D(2, -4, 0.0));
    qt3DWindow->camera()->setViewCenter(QVector3D(0, 0, 0));

    // For camera controls
    camController->setLinearSpeed(50.f);
    camController->setLookSpeed(180.f);
    camController->setCamera(qt3DWindow->camera());

    // Isosurface

    auto *entity = new Qt3DCore::QEntity(rootEntity);

    new Cylinder(entity, Qt::red,  {{-0.5f,0,0},{0.5f,0,0}}, 0.005f);
    new Cylinder(entity, Qt::green,{{0,-0.5f,0},{0,0.5f,0}}, 0.005f);
    new Cylinder(entity, Qt::blue, {{0,0,-0.5f},{0,0,0.5f}}, 0.005f);
    new Surface(entity, surfaceData, Qt::gray,0.25);

    widget->show();

    QApplication::exec();
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

    int argc = 0;
    char* argv[] = {};
    QApplication app(argc,argv);

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
    new Cylinder(entity, Qt::red,  {{-0.5f,0,0},{0.5f,0,0}}, 0.005f);
    new Cylinder(entity, Qt::green,{{0,-0.5f,0},{0,0.5f,0}}, 0.005f);
    new Cylinder(entity, Qt::blue, {{0,0,-0.5f},{0,0,0.5f}}, 0.005f);
    new Surface(entity, surfaceData, Qt::gray);

    qt3DWindow->setRootEntity(rootEntity);
    widget->show();

    QApplication::exec();
}