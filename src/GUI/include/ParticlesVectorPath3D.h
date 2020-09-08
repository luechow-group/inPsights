// Copyright (C) 2017-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_PARTICLECOLLECTIONPATH3D_H
#define INPSIGHTS_PARTICLECOLLECTIONPATH3D_H

#include <ParticlesVectorCollection.h>
#include <Qt3DCore/QEntity>

class ParticlesVectorPath3D{
public:

    ParticlesVectorPath3D(Qt3DCore::QEntity *root,
                             const ElectronsVectorCollection &electronsVectorCollection,
                             float radius = 0.01f);
};

#endif //INPSIGHTS_PARTICLECOLLECTIONPATH3D_H
