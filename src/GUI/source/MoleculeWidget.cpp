// Copyright (C) 2017-2019 Michael Heuer.
// Copyright (C) 2020 Leonard Reuter
// SPDX-License-Identifier: GPL-3.0-or-later

#include <QtWidgets>
#include <QScreen>
#include <MoleculeWidget.h>

#include <InPsightsWidget.h>
#include <Metrics.h>
#include <Line3D.h>
#include <Bond3D.h>
#include <ColorPalette.h>
#include <Surface.h>
#include <SurfaceDataGenerator.h>
#include <QuickHull.h>
#include <Conversion.h>
#include <SpinCorrelations3D.h>
#include <spdlog/spdlog.h>
#include <X3domConverter.h>

MoleculeWidget::MoleculeWidget(QWidget *parent)
        :
        QWidget(parent),
        compatabilityMode_(false),
        qt3DWindow_(new Qt3DExtras::Qt3DWindow()),
        root_(new Qt3DCore::QEntity()),
        moleculeEntity_(new Qt3DCore::QEntity(root_)),
        cameraController_(new Qt3DExtras::QOrbitCameraController(root_)),
        screenshotButton_(new QPushButton("Save image", this)),
        x3dExportButton_(new QPushButton("Save x3d", this)),
        sedsExportButton_(new QPushButton("Export SEDs", this)),
        resetCameraButton_(new QPushButton("Reset camera", this)),
        pan_(new QSpinBox(this)),
        tilt_(new QSpinBox(this)),
        roll_(new QSpinBox(this)),
        zoom_(new QSpinBox(this)),
        initPan_(0),
        initTilt_(0),
        initRoll_(0),
        initZoom_(100),
        defaultCameraRadius_(0.0f),
        fileInfoText_(new QLabel("Info text")),
        panTiltRollText_(new QLabel("Pan/Tilt/Roll")),
        zoomText_(new QLabel("Zoom")),
        atomsVector3D_(nullptr),
        cartesianAxes_(nullptr),
        light_(new Qt3DRender::QDirectionalLight(root_)){

    auto outerLayout = new QVBoxLayout(this);
    auto innerLayout = new QHBoxLayout();

    setLayout(outerLayout);
    outerLayout->addWidget(createWindowContainer(qt3DWindow_),1);
    outerLayout->addLayout(innerLayout,0);
    innerLayout->addWidget(screenshotButton_,4);
    innerLayout->addWidget(x3dExportButton_,4);
    innerLayout->addWidget(sedsExportButton_,4);
    innerLayout->addWidget(fileInfoText_, 10);
    fileInfoText_->setAlignment(Qt::AlignCenter);
    innerLayout->addWidget(resetCameraButton_,4);
    innerLayout->addWidget(zoomText_,1);
    innerLayout->addWidget(zoom_,2);
    innerLayout->addWidget(panTiltRollText_,2);
    innerLayout->addWidget(pan_,2);
    innerLayout->addWidget(tilt_,2);
    innerLayout->addWidget(roll_,2);

    qt3DWindow_->setRootEntity(root_);

    connect(screenshotButton_, &QPushButton::clicked, this, &MoleculeWidget::onScreenshot);
    connect(x3dExportButton_, &QPushButton::clicked, this, &MoleculeWidget::onX3dExport);
    connect(sedsExportButton_, &QPushButton::clicked, this, &MoleculeWidget::onSedsExport);
    connect(resetCameraButton_, &QPushButton::clicked, this, &MoleculeWidget::onResetCamera);

    connect(pan_, qOverload<int>(&QSpinBox::valueChanged),
            this, &MoleculeWidget::onCameraBoxesChanged);
    connect(tilt_, qOverload<int>(&QSpinBox::valueChanged),
            this, &MoleculeWidget::onCameraBoxesChanged);
    connect(roll_, qOverload<int>(&QSpinBox::valueChanged),
            this, &MoleculeWidget::onCameraBoxesChanged);
    connect(zoom_, qOverload<int>(&QSpinBox::valueChanged),
            this, &MoleculeWidget::onCameraBoxesChanged);

    // initialize light
    light_->setColor("white");
    light_->setIntensity(0.3);
    root_->addComponent(light_);

    setMouseTracking(true);
}

int MoleculeWidget::getAtomsNumber() {
    return sharedAtomsVector_->positionsVector().numberOfEntities();
}

void MoleculeWidget::onCameraBoxesChanged(int) {
    defaultCameraView();
    qt3DWindow_->camera()->panAboutViewCenter(float(pan_->value()));
    qt3DWindow_->camera()->tiltAboutViewCenter(float(tilt_->value()));
    qt3DWindow_->camera()->rollAboutViewCenter(float(roll_->value()));
    qt3DWindow_->camera()->viewSphere({0, 0, 0},
                                      100.0f / float(zoom_->value()) * defaultCameraRadius_);

    QVector3D v1 = -qt3DWindow_->camera()->position()/qt3DWindow_->camera()->position().length();
    QVector3D v2 = -qt3DWindow_->camera()->upVector();
    light_->setWorldDirection(v1 + v2 + QVector3D::crossProduct(v1,v2));
}

Qt3DCore::QEntity *MoleculeWidget::getMoleculeEntity() {
    return moleculeEntity_;
}

void MoleculeWidget::setupCameraBoxes(int pan, int tilt, int roll, int zoom) {
    pan_->setRange(-180,180);
    pan_->setWrapping(true);
    pan_->setSingleStep(5);
    pan_->setValue(pan);
    pan_->setSuffix(" °");
    pan_->setButtonSymbols(QAbstractSpinBox::PlusMinus);
    pan_->setAccelerated(true);

    tilt_->setRange(-180,180);
    tilt_->setWrapping(true);
    tilt_->setSingleStep(5);
    tilt_->setValue(tilt);
    tilt_->setSuffix(" °");
    tilt_->setButtonSymbols(QAbstractSpinBox::PlusMinus);
    tilt_->setAccelerated(true);

    roll_->setRange(-180,180);
    roll_->setWrapping(true);
    roll_->setSingleStep(5);
    roll_->setValue(roll);
    roll_->setSuffix(" °");
    roll_->setButtonSymbols(QAbstractSpinBox::PlusMinus);
    roll_->setAccelerated(true);

    // max value for zoom is 999 instead of INT_MAX, since it determines the box width
    zoom_->setRange(0,999);
    zoom_->setSingleStep(5);
    zoom_->setValue(zoom);
    zoom_->setSuffix(" %");
    zoom_->setButtonSymbols(QAbstractSpinBox::PlusMinus);
    zoom_->setAccelerated(true);
}

void MoleculeWidget::initialCameraSetup(int zoom, int pan, int tilt, int roll) {
    initZoom_ = zoom; initPan_ = pan; initTilt_ = tilt; initRoll_ = roll;

    calculateDefaultCameraRadius();

    cameraController_->setLinearSpeed(50.f);
    cameraController_->setLookSpeed(180.f);

    cameraController_->setCamera(qt3DWindow_->camera());

    resetCamera();
}

void MoleculeWidget::calculateDefaultCameraRadius() {
    float maxDistanceFromCenter = 0.0;
    for (Eigen::Index i = 0; i < atomsVector3D_->positionsVector().numberOfEntities(); ++i) {
        auto distanceFromCenter = static_cast<float>(sharedAtomsVector_->positionsVector()[i].norm())
                                  + GuiHelper::radiusFromType<Element>(sharedAtomsVector_->typesVector()[i]);
        if (distanceFromCenter > maxDistanceFromCenter)
            maxDistanceFromCenter = distanceFromCenter;
    }

    // add padding
    defaultCameraRadius_ += maxDistanceFromCenter + GuiHelper::radiusFromType<Element>(Element::H);
    spdlog::info("Determined camera radius with {} [a0]", defaultCameraRadius_);
}

void MoleculeWidget::resetCamera() {
    setupCameraBoxes(initPan_, initTilt_, initRoll_, initZoom_);
    onCameraBoxesChanged(0);
}

void MoleculeWidget::defaultCameraView() {
    qt3DWindow_->camera()->setViewCenter({0,0,0});
    qt3DWindow_->camera()->setUpVector({1,0,0});
    qt3DWindow_->camera()->setPosition({0,1,0});
    qt3DWindow_->camera()->viewSphere({0, 0, 0}, defaultCameraRadius_);
}

void MoleculeWidget::drawAxes(bool drawQ) {
    if (drawQ) {
        cartesianAxes_ = new CartesianAxes(moleculeEntity_);
    } else {
        cartesianAxes_->deleteLater();
    }
}

void MoleculeWidget::drawAtoms(bool drawQ) {
    if (drawQ) {
        atomsVector3D_ = new AtomsVector3D(moleculeEntity_, *sharedAtomsVector_);
    } else {
        atomsVector3D_->deleteLater();
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

void MoleculeWidget::addElectronsVector(const ElectronsVector &electronsVector, int clusterId, int structureId, bool coloredQ) {
    activeElectronsVectorsMap_[clusterId][structureId] = new ElectronsVector3D(moleculeEntity_, electronsVector, coloredQ);
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

void MoleculeWidget::addSeds(int clusterId, int structureId, const std::vector<ClusterData> &clusterData, double includedPercentage) {
    auto spins = clusterData[clusterId].exemplaricStructures_[structureId].typesVector();
    auto N = spins.numberOfEntities();

    std::vector<Surface*> seds(N);

    const auto& voxelData = clusterData[clusterId].voxelCubes_;

    for (long i = 0; i < spins.numberOfEntities(); ++i) {
        SurfaceDataGenerator surfaceDataGenerator(voxelData[i]);
        auto surfaceData = surfaceDataGenerator.computeSurfaceData(includedPercentage);
        seds[i] = new Surface(getMoleculeEntity(), surfaceData, GuiHelper::QColorFromType(spins[i]), 0.15);
    }
    activeSedsMap_[clusterId] = seds;
}

void MoleculeWidget::removeSeds(int clusterId) {
    for (auto &sed : activeSedsMap_[clusterId])
        sed->deleteLater();

    if (!activeSedsMap_[clusterId].empty())
        activeSedsMap_.erase(clusterId);
}


void MoleculeWidget::drawSpinCorrelations(const std::vector<ClusterData> &clusterData,
                                          double spinCorrelationThreshold, bool drawSameSpinCorrelationsQ) {
    for (auto &cluster : activeElectronsVectorsMap_)
        for (auto &structure : cluster.second) {
                new SpinCorrelations3D(structure.second,
                        clusterData[cluster.first].SeeStats_,
                        spinCorrelationThreshold,
                        drawSameSpinCorrelationsQ,
                        compatabilityMode_);
        }
}

void MoleculeWidget::deleteSpinCorrelations() {
    for (auto &cluster : activeElectronsVectorsMap_)
        for (auto &structure : cluster.second) {
            structure.second->deleteCorrelations();
        }
};

void MoleculeWidget::onAtomsChecked(int selectedParticle) {
    auto &particles = atomsVector3D_->particles3D_;

    for (int i = 0; i < static_cast<int>(particles.size()); ++i) {
        particles[i]->onSelected(i == selectedParticle);
    }
}

void MoleculeWidget::onAtomsHighlighted(int selectedParticle) {
    auto &particles = atomsVector3D_->particles3D_;

    for (int i = 0; i < static_cast<int>(particles.size()); ++i) {
        particles[i]->onHighlighted(i == selectedParticle);
    }
}

void MoleculeWidget::onElectronsChecked(int selectedParticle) {
    for (auto &cluster : activeElectronsVectorsMap_){
        auto &particles = cluster.second.begin()->second->particles3D_;
        for (int i = 0; i < static_cast<int>(particles.size()); ++i) {
            particles[i]->onSelected(i == selectedParticle);
        }
    }
}

void MoleculeWidget::onElectronsHighlighted(int selectedParticle) {
    for (auto &cluster : activeElectronsVectorsMap_){
        auto &particles = cluster.second.begin()->second->particles3D_;
        for (int i = 0; i < static_cast<int>(particles.size()); ++i) {
            particles[i]->onHighlighted(i == selectedParticle);
        }
    }
}

void MoleculeWidget::addMaximaHulls(int clusterId, const std::vector<ClusterData> &clusterData) {
    auto spins = clusterData[clusterId].representativeStructure().typesVector();
    auto N = spins.numberOfEntities();

    std::vector<Surface*> maximaHulls(N);

    quickhull::QuickHull<float> qh;

    for (long i = 0; i < spins.numberOfEntities(); ++i) {
        std::vector<Eigen::Vector3f> electronPointCloud;
        for(const auto & ev : clusterData[clusterId].exemplaricStructures_)
            electronPointCloud.emplace_back(ev[i].position().cast<float>());

        auto hull = qh.getConvexHull(electronPointCloud);
        auto vertices = hull.getVertices();
        auto triangles = hull.getTriangles();

        Conversion::calculateVertexNormals(vertices, triangles);
        SurfaceData surfaceData;
        surfaceData.triangles = triangles;
        surfaceData.vertices = vertices;

        maximaHulls[i] = new Surface(getMoleculeEntity(), surfaceData, ColorPalette::colorFunction(i), 0.5);
    }
    activeMaximaHullsMap_[clusterId] = maximaHulls;
}

void MoleculeWidget::removeMaximaHulls(int clusterId) {
    for (auto &hull : activeMaximaHullsMap_[clusterId])
        hull->deleteLater();

    if (!activeMaximaHullsMap_[clusterId].empty())
        activeMaximaHullsMap_.erase(clusterId);
}

void MoleculeWidget::onScreenshot(bool) {
    QScreen *screen = QGuiApplication::primaryScreen();
    if (const QWindow *window = windowHandle())
        screen = window->screen();
    if (!screen)
        return;
    std::string name = createFilenameFromActiveElectronvectors() + ".png";

    QFile file(name.c_str());
    file.open(QIODevice::WriteOnly);
    screen->grabWindow(qt3DWindow_->winId()).save(&file, "PNG");
    file.close();
}

void MoleculeWidget::onResetCamera(bool) {
    resetCamera();
}

std::string MoleculeWidget::createFilenameFromActiveElectronvectors() const {
    std::string name;

    // if  the inPsightsWidget is the parent
    auto inPsightsWidget = dynamic_cast<InPsightsWidget*>(parent());
    if(inPsightsWidget) {
        name = inPsightsWidget->filenameWithoutExtension().c_str();

        for (auto [key, submap] : activeElectronsVectorsMap_) {
            name += "-";
            name += std::to_string(key);
            name += "-[";
            if(inPsightsWidget->plotAllActiveQ()){
                name  += "all]";
            } else {
                for (auto[subkey, value] : submap) {
                    name += std::to_string(subkey);
                    name += ",";
                }
                name = name.substr(0, name.size() - 1);
                name += "]";
            }
        }
    } else {
        name = QDateTime::currentDateTime().toString(Qt::ISODate).toStdString();
    }
    return name;
}


void MoleculeWidget::onX3dExport(bool) {

    auto filename = createFilenameFromActiveElectronvectors() + ".html";

    X3domConverter x3Dconverter(filename, filename, light_->worldDirection());
    auto atoms3d = atomsVector3D_->particles3D_;

    for (const auto & a : atoms3d) {
        x3Dconverter.addSphere(a.operator*(),2);
    }

    //TODO Refactor
    double bondDrawingLimit = 1.40 * 1e-10 / AU::length;
    for (long i = 0; i < atomsVector3D_->numberOfEntities(); ++i) {
        for (long j = i+1; j < atomsVector3D_->numberOfEntities(); ++j) {

            auto atomDistance = Metrics::distance(
                    atomsVector3D_->operator[](i).position(),
                    atomsVector3D_->operator[](j).position());
            auto addedGuiRadii = (GuiHelper::radiusFromType(atomsVector3D_->operator[](i).type())
                                  + GuiHelper::radiusFromType(atomsVector3D_->operator[](j).type()));


            if (atomDistance - 0.5*addedGuiRadii < bondDrawingLimit) {
                Qt3DCore::QEntity root;
                Bond3D bond(&root, *atomsVector3D_->particles3D_[i], *atomsVector3D_->particles3D_[j]);
                x3Dconverter.addCylinder(*bond.srcCylinder_,2);
                x3Dconverter.addCylinder(*bond.destCylinder_,2);
            }
        }
    }

    for(const auto& evMap : activeElectronsVectorsMap_){
        Qt3DCore::QEntity root;

        for(const auto& ev : evMap.second) {
            for (const auto &e : ev.second->particles3D_) {
                x3Dconverter.addSphere(e.operator*(),1);
            }
        }
    }

    x3Dconverter.closeScene();
}

void MoleculeWidget::activateCompatabilityMode() {
    compatabilityMode_ = true;
}

void MoleculeWidget::onSedsExport(bool) {
    auto inPsightsWidget = dynamic_cast<InPsightsWidget*>(parent());
    if(inPsightsWidget) {
        for (auto &sed: activeSedsMap_) {
            auto clusterId = sed.first;
            int i = 0;
            for (auto &voxelCube: inPsightsWidget->clusterCollection_[clusterId].voxelCubes_) {
                std::string fileName = "sed-" + std::to_string(clusterId) + "-" + std::to_string(i) + ".txt";
                std::string comment = " cluster " + std::to_string(clusterId) + " SED " + std::to_string(i);
                voxelCube.exportMacmolplt(fileName, comment);
                i += 1;
            }
        }
        spdlog::info("Exported SEDs");
    }
}
