//
// Created by Michael Heuer on 12.11.17.
//

#ifndef AMOLQCGUI_PARTICLECOLLECTIONPATH3D_H
#define AMOLQCGUI_PARTICLECOLLECTIONPATH3D_H

#include "AtomCollection.h"
#include "ElectronCollections.h"

class ParticleCollectionPath3D{
public:

    ParticleCollectionPath3D(Qt3DCore::QEntity *root,
                             const ElectronCollections &electronCollections,
                             float radius = 0.01f);
};

#endif //AMOLQCGUI_PARTICLECOLLECTIONPATH3D_H
