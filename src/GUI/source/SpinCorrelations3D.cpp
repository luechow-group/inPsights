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

#include <SpinCorrelations3D.h>
#include <Line3D.h>
#include <Cylinder.h>

SpinCorrelations3D::SpinCorrelations3D(ElectronsVector3D *electronsVector3D,
                                       const TriangularMatrixStatistics &SeeStats,
                                       double spinCorrelationThreshold,
                                       bool drawSameSpinCorrelationsQ,
                                       bool compatabilityMode)
        :
        IConnection(electronsVector3D->correlations_),
        compatabilityMode_(compatabilityMode) {

    createConnections(*electronsVector3D, SeeStats, spinCorrelationThreshold, drawSameSpinCorrelationsQ);
}

void SpinCorrelations3D::createConnections(const ElectronsVector &electronsVector,
                                           const TriangularMatrixStatistics &SeeStats,
                                           double spinCorrelationThreshold,
                                           bool drawSameSpinCorrelationsQ) {

    auto electronRadius = float(GuiHelper::radiusFromType(Spin::alpha));

    // make list of electrons at same positions

    std::list<Eigen::Index> electronsAtSamePositions;
    if(compatabilityMode_) {
        for (auto i = 0; i < electronsVector.numberOfEntities() - 1; ++i) {
            for (auto j = i + 1; j < electronsVector.numberOfEntities(); ++j) {
                if (SpinPairClassification::atSamePositionQ({i, j}, electronsVector)) {
                    electronsAtSamePositions.push_back(i);
                    electronsAtSamePositions.push_back(j);
                }
            }
        }
    }

    // iterate over all pairs
    for (auto i = 0; i < electronsVector.numberOfEntities() - 1; ++i) {
        for (auto j = i + 1; j < electronsVector.numberOfEntities(); ++j) {

            // if i or j are not in the list of electrons at same positions
            if (std::find(electronsAtSamePositions.begin(), electronsAtSamePositions.end(), i) ==
                electronsAtSamePositions.end()
                && std::find(electronsAtSamePositions.begin(), electronsAtSamePositions.end(), j) ==
                   electronsAtSamePositions.end()) {

                const auto corr = float(SeeStats.mean()(i, j));
                const auto absCorr = std::abs(corr);

                if (absCorr >= spinCorrelationThreshold) {
                    Cylinder *c;

                    auto positionPair = GuiHelper::sphericalSurfacePositionPair(
                            electronsVector.positionsVector()[i], electronRadius,
                            electronsVector.positionsVector()[j], electronRadius);

                    if (corr < 0) {
                        c = new Cylinder(this, Qt::green, positionPair, electronRadius / 7.5f, absCorr);
                        c->material->setShininess(0);
                    } else if (drawSameSpinCorrelationsQ) {
                        c = new Cylinder(this, Qt::magenta, positionPair, electronRadius / 7.5f, absCorr);
                        c->material->setShininess(0);
                    }
                }
            }
        }
    }
}
