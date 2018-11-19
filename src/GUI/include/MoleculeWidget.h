//
// Created by Michael Heuer on 12.11.17.
//

#ifndef AMOLQCPP_MOLECULEWIDGET_H
#define AMOLQCPP_MOLECULEWIDGET_H

#include <QWidget>
#include <Qt3DCore>
#include <Qt3DExtras>
#include <QVBoxLayout>
#include <QLabel>

#include <AtomsVector3D.h>
#include <ElectronsVector3D.h>

class MoleculeWidget : public QWidget{
    Q_OBJECT
public:
    explicit MoleculeWidget(QWidget *parent = nullptr);
    Qt3DCore::QEntity* getRoot();

    void setElectrons(const ElectronsVector& electrons){
        //TODO delete old electronsVector3D_;
        electronsVector3D_ = new ElectronsVector3D(root_, electrons);
    }

    void setAtoms(const AtomsVector& atoms){
        //delete atomsVector3D_;
        atomsVector3D_ = new AtomsVector3D(root_,atoms);
    }

private:
    QVBoxLayout *layout_;
    Qt3DExtras::Qt3DWindow *qt3DWindow_;
    Qt3DCore::QEntity *root_;
    Qt3DExtras::QOrbitCameraController *cameraController_;
    QLabel* infoText_;


    AtomsVector3D* atomsVector3D_;
    ElectronsVector3D* electronsVector3D_;
};

#endif //AMOLQCPP_MOLECULEWIDGET_H
