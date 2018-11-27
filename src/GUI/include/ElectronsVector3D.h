//
// Created by Michael Heuer on 12.11.17.
//

#ifndef INPSIGHTS_ELECTRONCOLLECTION3D_H
#define INPSIGHTS_ELECTRONCOLLECTION3D_H

#include "Abstract3dObject.h"
#include <ParticlesVector.h>

class ElectronsVector3D : public Qt3DCore::QEntity {
public:
    ElectronsVector3D(Qt3DCore::QEntity *root,
            const ElectronsVector
            &electonsVector,
            bool showIndicesQ = false);

    ElectronsVector3D(Qt3DCore::QEntity *root,
            const AtomsVector& atomsVector,
            const ElectronsVector &electonsVector,
            bool showIndicesQ = false);

    void drawElectrons(bool showIndicesQ = false);

    void drawConnections(const AtomsVector& atomsVector,
            const ElectronsVector &electronsVector);

    ElectronsVector  electronsVector_;
};

#endif //INPSIGHTS_ELECTRONCOLLECTION3D_H
