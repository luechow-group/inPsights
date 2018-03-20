//
// Created by Michael Heuer on 12.11.17.
//

#ifndef AMOLQCPP_PARTICLECOLLECTIONPATH3D_H
#define AMOLQCPP_PARTICLECOLLECTIONPATH3D_H

#include "AtomsVector.h"
#include "ElectronsVectorCollection.h"

class ParticlesVectorPath3D{
public:

    ParticlesVectorPath3D(Qt3DCore::QEntity *root,
                             const ElectronsVectorCollection &electronsVectorCollection,
                             float radius = 0.01f);
};

#endif //AMOLQCPP_PARTICLECOLLECTIONPATH3D_H
