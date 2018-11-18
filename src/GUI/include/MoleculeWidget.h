//
// Created by Michael Heuer on 12.11.17.
//

#ifndef AMOLQCPP_MOLECULEWIDGET_H
#define AMOLQCPP_MOLECULEWIDGET_H

#include <Qt3DCore>
#include <Qt3DExtras>

class MoleculeWidget{
public:
    MoleculeWidget();
    Qt3DCore::QEntity* getRoot();
    QWidget* getWidget();

private:
    Qt3DCore::QEntity* root_;
    Qt3DExtras::Qt3DWindow* qt3DWindow_;
    QWidget* windowContainer_;
    Qt3DExtras::QOrbitCameraController* cameraController_;
};

#endif //AMOLQCPP_MOLECULEWIDGET_H
