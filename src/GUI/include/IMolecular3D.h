//
// Created by heuer on 30.11.18.
//

#ifndef INPSIGHTS_IMOLECULAR3D_H
#define INPSIGHTS_IMOLECULAR3D_H

#include <MoleculeWidget.h>

class IMolecular3D {
public:
    IMolecular3D(MoleculeWidget* moleculeWidgetParent);

    MoleculeWidget* moleculeWidgetParent_;
};

#endif //INPSIGHTS_IMOLECULAR3D_H
