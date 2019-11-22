/* Copyright (C) 2017-2019 Michael Heuer.
 *
 * This file is part of inPsights.
 * inPsights is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * inPsights is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with inPsights. If not, see <https://www.gnu.org/licenses/>.
 */

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
    void drawSpinConnections(bool drawQ = true);
    void drawSpinCorrelations(bool drawQ,
                              const std::vector<ClusterData> &clusterData,
                              double spinCorrelationThreshold, bool drawSameSpinCorrelationsQ);

    void initialCameraSetup(float distance = 8.0f,float pan = 0.0f, float tilt = 45.0f, float roll = 0.0f);
    void setupSpinBoxes(float pan, float tilt, float roll);
    void defaultCameraView();

    void setSharedAtomsVector(AtomsVector atomsVector);
    void addElectronsVector(const ElectronsVector& electronsVector, int clusterId = 0, int structureId = 0, bool coloredQ = false);
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

    void onCameraSpinBoxesChanged(int);
    void onScreenshot(bool);
    void onX3dExport(bool);

private:
    Qt3DExtras::Qt3DWindow *qt3DWindow_;
    Qt3DCore::QEntity *root_, *moleculeEntity_;
    Qt3DExtras::QOrbitCameraController *cameraController_;
    QPushButton *screenshotButton_, *x3dExportButton_;
    QSpinBox *pan_, *tilt_, *roll_;
    float defaultCameraDistance_;

public:
    QLabel* fileInfoText_, *panTiltRollText_;
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
