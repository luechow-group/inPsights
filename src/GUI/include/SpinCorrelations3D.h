/* Copyright (C) 2018-2019 Michael Heuer.
 *
 * This file is part of inPsights.
 * inPsights is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * inPsights is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with inPsights. If not, see <https://www.gnu.org/licenses/>.
 */

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
