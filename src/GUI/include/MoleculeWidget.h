//
// Created by Michael Heuer on 12.11.17.
//

#ifndef INPSIGHTS_MOLECULEWIDGET_H
#define INPSIGHTS_MOLECULEWIDGET_H

#include <QWidget>
#include <Qt3DCore>
#include <Qt3DExtras>
#include <QVBoxLayout>
#include <QLabel>
#include <QPushButton>
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
    void drawSpinConnections(bool drawQ = true);
    void drawSpinCorrelations(bool drawQ,
                              const std::vector<ClusterData> &clusterData,
                              double spinCorrelationThreshold);


    void setSharedAtomsVector(AtomsVector atomsVector);
    void addElectronsVector(const ElectronsVector& electronsVector, int clusterId = 0, int structureId = 0);
    void removeElectronsVector(int clusterId = 0, int structureId = 0);

    void addSeds(int clusterId, const std::vector<ClusterData> &clusterData, double includedPercentage);
    void removeSeds(int clusterId);

    void addMaximaHulls(int clusterId, const std::vector<ClusterData> &clusterData);
    void removeMaximaHulls(int clusterId);

public Q_SLOTS:

    void onAtomsChecked(std::vector<int>);
    void onElectronsChecked(std::vector<int>);
    void onAtomsHighlighted(std::vector<int>);
    void onElectronsHighlighted(std::vector<int>);

    void onScreenshot(bool);

private:
    Qt3DExtras::Qt3DWindow *qt3DWindow_;
    Qt3DCore::QEntity *root_, *moleculeEntity_;
    Qt3DExtras::QOrbitCameraController *cameraController_;
    QPushButton *screenshotButton_;
public:
    QLabel* infoText_;
private:
    std::shared_ptr<AtomsVector> sharedAtomsVector_;
    AtomsVector3D *atomsVector3D_;
    CartesianAxes *cartesianAxes_;

public:
    std::map<int, std::vector<Surface*>> activeSedsMap_, activeMaximaHullsMap_;
    std::map<int, std::map<int,ElectronsVector3D*>> activeElectronsVectorsMap_;

};

#endif //INPSIGHTS_MOLECULEWIDGET_H
