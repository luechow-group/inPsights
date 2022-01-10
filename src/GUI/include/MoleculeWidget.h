// Copyright (C) 2017-2019 Michael Heuer.
// Copyright (C) 2020 Leonard Reuter
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_MOLECULEWIDGET_H
#define INPSIGHTS_MOLECULEWIDGET_H

#include <QWidget>
#include <Qt3DCore>
#include <Qt3DExtras>
#include <QVBoxLayout>
#include <QLabel>
#include <QPushButton>
#include <QSpinBox>
#include <ParticlesVector3D.h>
#include <Statistics.h>
#include <ClusterData.h>
#include <SurfaceData.h>
#include <Surface.h>
#include <Polyline.h>
#include <CartesianAxes.h>

class MoleculeWidget : public QWidget{
    Q_OBJECT
public:
    explicit MoleculeWidget(QWidget *parent = nullptr);

    Qt3DCore::QEntity* getMoleculeEntity();

    //TODO make base MoleculeWidget and InPsightsMoleculeWidget child

    void drawAxes(bool drawQ = true);
    void drawAtoms(bool drawQ = true);
    void drawBonds(bool drawQ = true);
    void drawSpinCorrelations(const std::vector<ClusterData> &clusterData,
                              double spinCorrelationThreshold, bool drawSameSpinCorrelationsQ);
    void deleteSpinCorrelations();

    void initialCameraSetup(int zoom, int pan, int tilt, int roll);
    void setupCameraBoxes(int pan, int tilt, int roll, int zoom);
    void defaultCameraView();
    void resetCamera();

    void calculateDefaultCameraRadius();

    void setSharedAtomsVector(AtomsVector atomsVector);
    void addElectronsVector(const ElectronsVector& electronsVector, int clusterId = 0, int structureId = 0, bool coloredQ = false);
    void removeElectronsVector(int clusterId = 0, int structureId = 0);

    void addSeds(int clusterId, int structureId, const std::vector<ClusterData> &clusterData, double includedPercentage);
    void removeSeds(int clusterId);

    void addMaximaHulls(int clusterId, const std::vector<ClusterData> &clusterData);
    void removeMaximaHulls(int clusterId);

    int getAtomsNumber();

public slots:
    void activateCompatabilityMode();
    void onAtomsChecked(int);
    void onElectronsChecked(int);
    void onAtomsHighlighted(int);
    void onElectronsHighlighted(int);

    void onCameraBoxesChanged(int);
    void onScreenshot(bool);
    void onX3dExport(bool);
    void onResetCamera(bool);

private:
    bool compatabilityMode_;
    Qt3DExtras::Qt3DWindow *qt3DWindow_;
    Qt3DCore::QEntity *root_, *moleculeEntity_;
    Qt3DExtras::QOrbitCameraController *cameraController_;
    QPushButton *screenshotButton_, *x3dExportButton_, *resetCameraButton_;
    QSpinBox *pan_, *tilt_, *roll_, *zoom_;
    Qt3DRender::QDirectionalLight *light_;
    int initPan_, initTilt_, initRoll_,initZoom_;
    float defaultCameraRadius_;
public:
    QLabel* fileInfoText_, *panTiltRollText_, *zoomText_;
private:
    std::shared_ptr<AtomsVector> sharedAtomsVector_;
    AtomsVector3D *atomsVector3D_;
    CartesianAxes *cartesianAxes_;
public:
    std::map<int, std::vector<Surface*>> activeSedsMap_, activeMaximaHullsMap_;
    std::map<int, std::map<int,ElectronsVector3D*>> activeElectronsVectorsMap_;

    std::string createFilenameFromActiveElectronvectors() const;


};

#endif //INPSIGHTS_MOLECULEWIDGET_H
