//
// Created by Michael Heuer on 12.11.17.
//

#include <Qt3DRender>
#include <Qt3DExtras>
#include <QtWidgets/QApplication>
#include <QtWidgets>

#include "MoleculeWidget.h"


Qt3DCore::QEntity* MoleculeWidget::createMoleculeWidget() {

    auto *root = new Qt3DCore::QEntity();

    auto *view = new Qt3DExtras::Qt3DWindow();
    view->defaultFrameGraph()->setClearColor(Qt::white);

    QWidget *moleculeView = QWidget::createWindowContainer(view);

    // camera
    view->camera()->lens()->setPerspectiveProjection(45.0f, 16.0f / 9.0f, 0.1f, 100.0f);
    view->camera()->setPosition(QVector3D(2.5, -8, 0.0));
    view->camera()->setViewCenter(QVector3D(0, 0, 0));

    // manipulator
    auto *manipulator = new Qt3DExtras::QOrbitCameraController(root);
    manipulator->setLinearSpeed(50.f);
    manipulator->setLookSpeed(180.f);
    manipulator->setCamera(view->camera());

    view->setRootEntity(root);


    moleculeView->resize(800, 800);
    moleculeView->show();


    return root;
}