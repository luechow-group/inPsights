// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_SPINCORRELATIONS3D_H
#define INPSIGHTS_SPINCORRELATIONS3D_H

#include "ParticlesVector3D.h"
#include "SpinPairClassification.h"
#include "IConnection.h"
#include <Statistics.h>


class SpinCorrelations3D : public IConnection {
public:
    SpinCorrelations3D(ElectronsVector3D *electronsVector3D,
                       const TriangularMatrixStatistics& SeeStats,
                       double spinCorrelationThreshold,
                       bool drawSameSpinCorrelationsQ = false,
                       bool compatabilityMode = false);

    void createConnections(const ElectronsVector &electronsVector,
                           const TriangularMatrixStatistics &SeeStats,
                           double spinCorrelationThreshold,
                           bool drawSameSpinCorrelationsQ);
private:
    bool compatabilityMode_;
};

#endif //INPSIGHTS_SPINCORRELATIONS3D_H
