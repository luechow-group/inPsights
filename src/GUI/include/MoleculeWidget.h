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
#include <ParticlesVector3D.h>
#include <Statistics.h>

class MoleculeWidget : public QWidget{
    Q_OBJECT
public:
    explicit MoleculeWidget(QWidget *parent = nullptr);

    Qt3DCore::QEntity* getRoot();

        //void mouseMoveEvent(QMouseEvent* event) override;

    void drawAtoms(bool drawQ = true);

    void drawSpinConnections(bool drawQ);

    void drawSpinCorrelations(bool drawQ,
            const Statistics::RunningStatistics<Eigen::MatrixXd, unsigned, true> &SeeStats,
                              int spinCorrelationThreshold);

    void drawBonds(bool drawQ);

    void addElectronsVector(const ElectronsVector& electronsVector, int id = 0);

    void removeElectronsVector(int id = 0);

    void setSharedAtomsVector(AtomsVector atomsVector);

private:
    QVBoxLayout *layout_;
    Qt3DExtras::Qt3DWindow *qt3DWindow_;
    Qt3DCore::QEntity *root_, *moleculeEntity_;
    Qt3DExtras::QOrbitCameraController *cameraController_;
public:
    QLabel* infoText_;
private:
    std::shared_ptr<AtomsVector> sharedAtomsVector_;
    AtomsVector3D* atomsVector3D_;
    std::map<int, ElectronsVector3D*> activeElectronsVectorsMap_;
};

#endif //INPSIGHTS_MOLECULEWIDGET_H
