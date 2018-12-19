//
// Created by Michael Heuer on 04.12.18.
//

#ifndef INPSIGHTS_SAMESPINCONNECTIONS3D_H
#define INPSIGHTS_SAMESPINCONNECTIONS3D_H

#include "IConnection.h"
#include "ParticlesVector3D.h"

class SameSpinConnections3D : public IConnection{
public:
    explicit SameSpinConnections3D(ElectronsVector3D *electronsVector3D,
            double maxDistance = 1.6, double identicalThreshold = 0.01);

    void createConnections(ElectronsVector3D *electronsVector3D,
                           double maxDistance, double identicalThreshold);
};

#endif //INPSIGHTS_SAMESPINCONNECTIONS3D_H
