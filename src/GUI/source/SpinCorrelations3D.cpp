// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <SpinCorrelations3D.h>
#include <Line3D.h>
#include <Cylinder.h>

SpinCorrelations3D::SpinCorrelations3D(ElectronsVector3D *electronsVector3D,
                                       const TriangularMatrixStatistics &SeeStats,
                                       double spinCorrelationThreshold,
                                       bool drawSameSpinCorrelationsQ,
                                       bool compatibilityMode)
        :
        IConnection(electronsVector3D->correlations_),
        compatibilityMode_(compatibilityMode) {

    createConnections(*electronsVector3D, SeeStats, spinCorrelationThreshold, drawSameSpinCorrelationsQ);
}

void SpinCorrelations3D::createConnections(const ElectronsVector &electronsVector,
                                           const TriangularMatrixStatistics &SeeStats,
                                           double spinCorrelationThreshold,
                                           bool drawSameSpinCorrelationsQ) {

    auto electronRadius = float(GuiHelper::radiusFromType(Spin::alpha));

    // make list of electrons at same positions

    std::list<Eigen::Index> electronsAtSamePositions;
    if(compatibilityMode_) {
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

                    Qt::GlobalColor color;
                    if (corr < 0) {
                        if (drawSameSpinCorrelationsQ) {
                            color = Qt::green;
                        }
                        else {
                            color = Qt::black;
                        }
                    } else {
                        if (drawSameSpinCorrelationsQ) {
                            color = Qt::magenta;
                        }
                        else {
                            continue;
                        }
                    }
                    if (positionPair.first != positionPair.second) {
                        c = new Cylinder(this, color, positionPair, electronRadius / 5.0f, absCorr);
                        c->material->setShininess(0);
                    }
                }
            }
        }
    }
}
