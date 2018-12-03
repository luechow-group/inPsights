//
// Created by Michael Heuer on 12.11.17.
//

#include <Qt3DRender>
#include <QtWidgets>
#include <MoleculeWidget.h>

#include <Metrics.h>
#include <Line3D.h>
#include <Bond3D.h>
#include <Cylinder.h>
#include <Particle3D.h>
#include <AtomsVectorLinkedElectronsVector.h>

#include <SpinCorrelations3D.h>

MoleculeWidget::MoleculeWidget(QWidget *parent)
    :
    QWidget(parent),
    layout_(new QVBoxLayout(this)),
    qt3DWindow_(new Qt3DExtras::Qt3DWindow()),
    root_(new Qt3DCore::QEntity()),
    moleculeEntity_(new Qt3DCore::QEntity(root_)),
    cameraController_(new Qt3DExtras::QOrbitCameraController(root_)),
    infoText_(new QLabel("Info text")),
    atomsVector3D_(nullptr)
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

    //setToolTip(QString("Hey there"));
    //setMouseTracking(true);
}

Qt3DCore::QEntity* MoleculeWidget::getRoot() {
    return root_;
}

/*void MoleculeWidget::mouseMoveEvent(QMouseEvent *event) {
    std::cout << "MOUSE MOVEMENT" << std::endl;

    QWidget::mouseMoveEvent(event);  // Or whatever the base class is.
}*/

void MoleculeWidget::drawAtoms(bool drawQ) {
    if(drawQ)
        atomsVector3D_ = new AtomsVector3D(moleculeEntity_, *sharedAtomsVector_);
    else
        atomsVector3D_->deleteLater();
}

void MoleculeWidget::drawBonds(bool drawQ) {
    if(drawQ)
        atomsVector3D_->drawConnections();
    else
        atomsVector3D_->deleteConnections();
}

void MoleculeWidget::drawSpinConnections(bool drawQ) {
    if(drawQ)
        for(auto& mapItem : activeElectronsVectorsMap_)
            mapItem.second->drawConnections();
    else
        for(auto& mapItem : activeElectronsVectorsMap_)
            mapItem.second->deleteConnections();
}

void MoleculeWidget::addElectronsVector(const ElectronsVector &electronsVector, int id) {
    activeElectronsVectorsMap_.emplace(id,new ElectronsVector3D(moleculeEntity_, electronsVector));
}

void MoleculeWidget::removeElectronsVector(int id) {
    activeElectronsVectorsMap_[id]->deleteLater();
    activeElectronsVectorsMap_.erase(id);
}

void MoleculeWidget::setSharedAtomsVector(AtomsVector atomsVector) {
    sharedAtomsVector_ = std::make_shared<AtomsVector >(std::move(atomsVector));
}

void MoleculeWidget::drawSpinCorrelations(bool drawQ,
                                          const Statistics::RunningStatistics<Eigen::MatrixXd, unsigned, true> &SeeStats,
                                          int spinCorrelationThreshold) {
  for(auto &mapItem : activeElectronsVectorsMap_) {
      if(drawQ) {
          new SpinCorrelations3D(mapItem.second, SeeStats, spinCorrelationThreshold);
      } else {
          std::cout << "off" << std::endl;
          mapItem.second->correlations_->deleteLater();
          mapItem.second->correlations_ = new Qt3DCore::QEntity(root_);
          //mapItem.second->deleteCorrelations();
          //activeElectronsVectorsMap_[mapItem.first]->deleteCorrelations();
      }
      //auto& connections =  electronsVector3D->iConnections_;
      /*
      if (drawQ) {
          connections.emplace_back(
                  new SpinCorrelations3D(electronsVector3D, SeeStats, spinCorrelationThreshold)
                  );
      } else {
          for(auto it = connections.begin(); it !=  connections.end(); it++) {
              // check if SpinCorrelations3D
              auto castedIConnection = dynamic_cast<SpinCorrelations3D*>(*it);
              if (castedIConnection) {
                  castedIConnection->deleteLater();
                  connections->erase(it);
              }
          }
      }*/
  }
}
