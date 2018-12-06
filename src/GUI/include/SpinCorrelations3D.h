#include <utility>

//
// Created by heuer on 03.12.18.
//

#ifndef INPSIGHTS_SPINCORRELATIONS3D_H
#define INPSIGHTS_SPINCORRELATIONS3D_H

#include <SpinConnections3D.h>
#include <Statistics.h>

class SpinCorrelations3D : public SpinConnections3D {
public:
    SpinCorrelations3D(ElectronsVector3D *electronsVector3D,
                       IntraParticlesStatistics SeeStats,
                       double spinCorrelationThreshold);

    //TODO add createConnection() as in SpinConnections3D, add update() method
    void drawSpinCorrelations(ElectronsVector3D *electronsVector3D, double spinCorrelationThreshold);

private:
    IntraParticlesStatistics SeeStats_;
};

#endif //INPSIGHTS_SPINCORRELATIONS3D_H
