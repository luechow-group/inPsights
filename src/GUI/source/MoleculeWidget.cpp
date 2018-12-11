//
// Created by Michael Heuer on 12.11.17.
//

#include <Qt3DRender>
#include <QtWidgets>
#include <QPixmap>
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
        qt3DWindow_(new Qt3DExtras::Qt3DWindow()),
        root_(new Qt3DCore::QEntity()),
        moleculeEntity_(new Qt3DCore::QEntity(root_)),
        cameraController_(new Qt3DExtras::QOrbitCameraController(root_)),
        screenshotButton_(new QPushButton("Save image", this)),
        infoText_(new QLabel("Info text")),
        atomsVector3D_(nullptr) {

    auto outerLayout = new QVBoxLayout(this);
    auto innerLayout = new QHBoxLayout();

    setLayout(outerLayout);
    outerLayout->addWidget(createWindowContainer(qt3DWindow_));
    outerLayout->addLayout(innerLayout);
    innerLayout->addWidget(screenshotButton_);
    innerLayout->addWidget(infoText_);

    infoText_->setFixedHeight(18);

    qt3DWindow_->setRootEntity(root_);

    qt3DWindow_->camera()->lens()->setPerspectiveProjection(45.0f, 16.0f / 9.0f, 0.1f, 100.0f);
    qt3DWindow_->camera()->setPosition(QVector3D(2.5, -8, 0.0));
    qt3DWindow_->camera()->setViewCenter(QVector3D(0, 0, 0));

    cameraController_->setLinearSpeed(50.f);
    cameraController_->setLookSpeed(180.f);
    cameraController_->setCamera(qt3DWindow_->camera());

    connect(screenshotButton_, &QPushButton::clicked, this, &MoleculeWidget::onScreenshot);

    setMouseTracking(true);
}

Qt3DCore::QEntity *MoleculeWidget::getMoleculeEntity() {
    return moleculeEntity_;
}

void MoleculeWidget::drawAtoms(bool drawQ) {
    if (drawQ) {
        atomsVector3D_ = new AtomsVector3D(moleculeEntity_, *sharedAtomsVector_);
    } else {
        atomsVector3D_->deleteConnections();
        atomsVector3D_->deleteLater();
        delete atomsVector3D_;
    }
}

void MoleculeWidget::drawBonds(bool drawQ) {
    if (atomsVector3D_) {
        if (drawQ)
            atomsVector3D_->drawConnections();
        else
            atomsVector3D_->deleteConnections();
    }
}

void MoleculeWidget::drawSpinConnections(bool drawQ) {
    if (drawQ)
        for (auto &cluster : activeElectronsVectorsMap_)
            for (auto &structure : cluster.second)
                structure.second->drawConnections();
    else
        for (auto &cluster : activeElectronsVectorsMap_)
            for (auto &structure : cluster.second)
                structure.second->deleteConnections();
}

void MoleculeWidget::addElectronsVector(const ElectronsVector &electronsVector, int clusterId, int structureId) {
    activeElectronsVectorsMap_[clusterId][structureId] = new ElectronsVector3D(moleculeEntity_, electronsVector);
}

void MoleculeWidget::removeElectronsVector(int clusterId, int structureId) {
    activeElectronsVectorsMap_[clusterId][structureId]->deleteLater();
    activeElectronsVectorsMap_[clusterId].erase(structureId);

    if(activeElectronsVectorsMap_[clusterId].empty()) {
        activeElectronsVectorsMap_.erase(clusterId);
    }
}

void MoleculeWidget::setSharedAtomsVector(AtomsVector atomsVector) {
    sharedAtomsVector_ = std::make_shared<AtomsVector>(std::move(atomsVector));
}

void MoleculeWidget::drawSpinCorrelations(bool drawQ,
                                          const std::vector<ClusterData> &clusterData,
                                          double spinCorrelationThreshold) {
    //TODO SPLIT INTO DRAW AND DELETE METHODS
    for (auto &cluster : activeElectronsVectorsMap_)
        for (auto &structure : cluster.second) {
            if (drawQ) {
                new SpinCorrelations3D(structure.second, clusterData[cluster.first].SeeStats_,
                                       spinCorrelationThreshold);
            } else {
                structure.second->correlations_->deleteLater();
                structure.second->correlations_ = new Qt3DCore::QEntity(structure.second);
            }
        }
}

void MoleculeWidget::onAtomsChecked(std::vector<int> selectedParticles) {
    auto &particles = atomsVector3D_->particles3D_;

    for (int i = 0; i < static_cast<int>(particles.size()); ++i) {
        auto foundQ = std::find(selectedParticles.begin(), selectedParticles.end(), i) != selectedParticles.end();
        particles[i]->onSelected(foundQ);
    }
}

void MoleculeWidget::onAtomsHighlighted(std::vector<int> selectedParticles) {
    auto &particles = atomsVector3D_->particles3D_;

    for (int i = 0; i < static_cast<int>(particles.size()); ++i) {
        auto foundQ = std::find(selectedParticles.begin(), selectedParticles.end(), i) != selectedParticles.end();
        particles[i]->onHighlighted(foundQ); //TODO add highlight, select and normal function
    }
}


void MoleculeWidget::onElectronsChecked(std::vector<int> selectedParticles) {
    if(activeElectronsVectorsMap_.size() != 1)
        std::cout << "Make sure only one Electronsvector is checked" << std::endl;
    else {
        auto &particles = activeElectronsVectorsMap_.begin()->second.begin()->second->particles3D_;

        for (int i = 0; i < static_cast<int>(particles.size()); ++i) {
            auto foundQ = std::find(selectedParticles.begin(), selectedParticles.end(), i) != selectedParticles.end();
            particles[i]->onSelected(foundQ);
        }
    }
}

void MoleculeWidget::onElectronsHighlighted(std::vector<int> selectedParticles) {
    if(activeElectronsVectorsMap_.size() != 1)
        std::cout << "Make sure only one Electronsvector is checked" << std::endl;
    else {
        auto &particles = activeElectronsVectorsMap_.begin()->second.begin()->second->particles3D_;

        for (int i = 0; i < static_cast<int>(particles.size()); ++i) {
            auto foundQ = std::find(selectedParticles.begin(), selectedParticles.end(), i) != selectedParticles.end();
            particles[i]->onHighlighted(foundQ);
        }
    }
}

void MoleculeWidget::onScreenshot(bool) {
    QScreen *screen = QGuiApplication::primaryScreen();
    if (const QWindow *window = windowHandle())
        screen = window->screen();
    if (!screen)
        return;

    QFile file(QDateTime::currentDateTime().toString(Qt::ISODate) + QString(".png"));
    file.open(QIODevice::WriteOnly);
    screen->grabWindow(qt3DWindow_->winId()).save(&file, "PNG");
    file.close();
}
