//
// Created by Michael Heuer on 12.11.17.
//

#ifndef AMOLQCPP_ELECTRONCOLLECTION3D_H
#define AMOLQCPP_ELECTRONCOLLECTION3D_H

#include "Abstract3dObject.h"
#include "ElectronsVector.h"
#include "AtomsVector.h"

class ElectronsVector3D{
public:
    ElectronsVector3D(Qt3DCore::QEntity *root, const ElectronsVector &electonsVector,
                             bool showIndicesQ = false);

    ElectronsVector3D(Qt3DCore::QEntity *root,
                      const AtomsVector& atomsVector,
                      const ElectronsVector &electonsVector,
                      bool showIndicesQ = false);

    void drawElectrons(Qt3DCore::QEntity *root, const ElectronsVector &electronsVector,bool showIndicesQ = false);

    void drawConnections(Qt3DCore::QEntity *root,
                         const AtomsVector& atomsVector,
                         const ElectronsVector &electronsVector);

};

#endif //AMOLQCPP_ELECTRONCOLLECTION3D_H
