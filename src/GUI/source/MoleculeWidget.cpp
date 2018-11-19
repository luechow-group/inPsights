//
// Created by Michael Heuer on 12.11.17.
//

#include <Qt3DRender>
#include <QtWidgets>

#include <MoleculeWidget.h>


MoleculeWidget::MoleculeWidget(QWidget *parent)
    :
    QWidget(parent),
    layout_(new QVBoxLayout(this)),
    qt3DWindow_(new Qt3DExtras::Qt3DWindow()),
    root_(new Qt3DCore::QEntity()),
    moleculeEntity_(new Qt3DCore::QEntity(root_)),
    cameraController_(new Qt3DExtras::QOrbitCameraController(root_)),
    infoText_(new QLabel("Info text")),
    atomsVector3D_(nullptr),
    electronsVector3D_(nullptr)
{
    setLayout(layout_);
    layout_->addWidget(createWindowContainer(qt3DWindow_));
    layout_->addWidget(infoText_);

    infoText_->setFixedHeight(14);

    qt3DWindow_->setRootEntity(root_);

    qt3DWindow_->camera()->lens()->setPerspectiveProjection(45.0f, 16.0f / 9.0f, 0.1f, 100.0f);
    qt3DWindow_->camera()->setPosition(QVector3D(2.5, -8, 0.0));
    qt3DWindow_->camera()->setViewCenter(QVector3D(0, 0, 0));

    cameraController_->setLinearSpeed(50.f);
    cameraController_->setLookSpeed(180.f);
    cameraController_->setCamera(qt3DWindow_->camera());
}

Qt3DCore::QEntity* MoleculeWidget::getRoot() {
    return root_;
}
