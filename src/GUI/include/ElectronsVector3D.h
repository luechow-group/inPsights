//
// Created by Michael Heuer on 12.11.17.
//

#ifndef AMOLQCPP_ELECTRONCOLLECTION3D_H
#define AMOLQCPP_ELECTRONCOLLECTION3D_H

#include "Abstract3dObject.h"
#include "ElectronsVector.h"

class ElectronsVector3D{
public:
    ElectronsVector3D(Qt3DCore::QEntity *root, const ElectronsVector &electonCollection,
                             bool showIndicesQ = false);
};

#endif //AMOLQCPP_ELECTRONCOLLECTION3D_H
