//
// Created by heuer on 03.12.18.
//

#ifndef INPSIGHTS_SPINCORRELATIONS3D_H
#define INPSIGHTS_SPINCORRELATIONS3D_H

#include "ParticlesVector3D.h"
#include "SpinClassification.h"
#include "IConnection.h"
#include <Statistics.h>


class SpinCorrelations3D : public IConnection {
public:
    SpinCorrelations3D(ElectronsVector3D *electronsVector3D,
                       const IntraParticlesStatistics& SeeStats,
                       double spinCorrelationThreshold);

    void createConnections(ElectronsVector3D *electronsVector3D,
                           const IntraParticlesStatistics &SeeStats,
                           double spinCorrelationThreshold);
};

#endif //INPSIGHTS_SPINCORRELATIONS3D_H
