//
// Created by Michael Heuer on 12.11.17.
//

#ifndef INPSIGHTS_PARTICLECOLLECTIONPATH3D_H
#define INPSIGHTS_PARTICLECOLLECTIONPATH3D_H

#include <ParticlesVectorCollection.h>

class ParticlesVectorPath3D{
public:

    ParticlesVectorPath3D(Qt3DCore::QEntity *root,
                             const ElectronsVectorCollection &electronsVectorCollection,
                             float radius = 0.01f);
};

#endif //INPSIGHTS_PARTICLECOLLECTIONPATH3D_H
