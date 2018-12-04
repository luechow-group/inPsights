//
// Created by Michael Heuer on 04.12.18.
//

#ifndef INPSIGHTS_SAMESPINCONNECTIONS3D_H
#define INPSIGHTS_SAMESPINCONNECTIONS3D_H

#include <SpinConnections3D.h>

class SameSpinConnections3D : public SpinConnections3D{
public:
    explicit SameSpinConnections3D(ElectronsVector3D *electronsVector3D,
            double identicalThreshold = 0.01, double maxDistance = 1.6);

    void createConnections(ElectronsVector3D *electronsVector3D);

};

#endif //INPSIGHTS_SAMESPINCONNECTIONS3D_H
