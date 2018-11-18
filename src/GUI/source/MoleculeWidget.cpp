//
// Created by Michael Heuer on 12.11.17.
//

#include <QWidget>
#include <MoleculeWidget.h>

MoleculeWidget::MoleculeWidget()
    :
    root_(new Qt3DCore::QEntity()),
    qt3DWindow_(new Qt3DExtras::Qt3DWindow()),
    windowContainer_(QWidget::createWindowContainer(qt3DWindow_)),
    cameraController_(new Qt3DExtras::QOrbitCameraController(root_))
{
    qt3DWindow_->setRootEntity(root_);

    // camera
    qt3DWindow_->camera()->lens()->setPerspectiveProjection(45.0f, 16.0f / 9.0f, 0.1f, 100.0f);
    qt3DWindow_->camera()->setPosition(QVector3D(2.5, -8, 0.0));
    qt3DWindow_->camera()->setViewCenter(QVector3D(0, 0, 0));

    cameraController_->setLinearSpeed(50.f);
    cameraController_->setLookSpeed(180.f);
    cameraController_->setCamera(qt3DWindow_->camera());

    //windowContainer_->resize(800, 800);
    //windowContainer_->show();
}

Qt3DCore::QEntity* MoleculeWidget::getRoot() {
    return root_;
}

QWidget *MoleculeWidget::getWidget() {
    return windowContainer_;
}
