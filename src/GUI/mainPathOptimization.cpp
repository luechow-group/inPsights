//
// Created by Michael Heuer on 07.11.17.
//

#include <iostream>

#include "OptimizationPathFileImporter.h"
#include "ElementInfo.h"

//#include <QGuiApplication>
#include <Qt3DCore>
#include <Qt3DRender>
#include <Qt3DExtras>

#include <QtWidgets/QApplication>
#include <QtWidgets/QMainWindow>
#include <QtCharts/QChartView>
#include <QtCharts/QLineSeries>
#include <QtCharts/QValueAxis>
#include <QtWidgets/QHBoxLayout>
#include <Qt3DExtras>

#include "MolecularGeometry3D.h"

#include "ElementInfo.h"
#include "ElementType.h"
#include "Sphere.h"
#include "Cylinder.h"
#include "DividedCylinder.h"
#include "Atom3D.h"
#include "Bond3D.h"
#include "Electron3D.h"
#include "Polyline.h"

#include <numeric>

int main(int argc, char *argv[]) {
    std::cout << "Hello" << std::endl;


    OptimizationPathFileImporter importer("Diborane-Paths.300",1);

    auto ecs = importer.getPath(1);

    std::cout << ecs.getSpinTypeCollection().spinTypesAsEigenVector().transpose() << std::endl;
    for (int i = 0; i < ecs.getNumberOfParticleCollections(); ++i) {
        std::cout << ecs[i].positionsAsEigenVector().transpose() << std::endl;
    }

    std::string str = "Hello World!";
    std::vector<std::string> vec(10,str);
    std::string a = std::accumulate(vec.begin(), vec.end(), std::string(""));
    std::cout << a << std::endl;

/*

    //MolecularGeometry3D molecularGeometry3D (root, waveFunctionParser.getAtomCollection());

    Qt3DCore::QEntity *root = new Qt3DCore::QEntity();

    Qt3DExtras::Qt3DWindow *view = new Qt3DExtras::Qt3DWindow();
    view->defaultFrameGraph()->setClearColor(Qt::white);

    QWidget *moleculeView = QWidget::createWindowContainer(view);

    Qt3DCore::QEntity *scene = root;
    // camera
    Qt3DRender::QCamera *camera = view->camera();
    camera->lens()->setPerspectiveProjection(45.0f, 16.0f / 9.0f, 0.1f, 100.0f);
    //camera->setPosition(QVector3D(4.25, 4.25, 0.0)); // ethane
    camera->setPosition(QVector3D(2.5, -5.00, 0.0)); // ethane
    camera->setViewCenter(QVector3D(0, 0, 0));
    //camera->setFarPlane(1000);
    //camera->setNearPlane(0.1);
    //camera->tilt(-5.0f);

    // manipulator
    //Qt3DExtras::QFirstPersonCameraController* manipulator = new Qt3DExtras::QFirstPersonCameraController(scene);
    Qt3DExtras::QOrbitCameraController *manipulator = new Qt3DExtras::QOrbitCameraController(scene);
    manipulator->setLinearSpeed(50.f);
    manipulator->setLookSpeed(180.f);
    manipulator->setCamera(camera);

    view->setRootEntity(scene);*/

};