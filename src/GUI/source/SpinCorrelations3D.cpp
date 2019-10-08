//
// Created by heuer on 03.12.18.
//

#include <SpinCorrelations3D.h>
#include <Line3D.h>
#include <Cylinder.h>

SpinCorrelations3D::SpinCorrelations3D(ElectronsVector3D *electronsVector3D,
                                       const TriangularMatrixStatistics& SeeStats,
                                       double spinCorrelationThreshold)
        :
        IConnection(electronsVector3D->correlations_) {

    createConnections(*electronsVector3D, SeeStats, spinCorrelationThreshold);
}

void SpinCorrelations3D::createConnections(const ElectronsVector &electronsVector,
                                           const TriangularMatrixStatistics &SeeStats,
                                           double spinCorrelationThreshold)  {

    auto pairTypes = SpinPairClassification::classify(electronsVector);
    auto electronRadius = float(GuiHelper::radiusFromType(Spin::alpha));

    for (auto &idxPair : pairTypes) {
        auto i = idxPair.first.first;
        auto j = idxPair.first.second;

        if (idxPair.second == SpinPairClassification::PairType::closeBy
        && !SpinPairClassification::isAtSamePositionQ(pairTypes,i)
        && !SpinPairClassification::isAtSamePositionQ(pairTypes,j)
        ) {
            auto corr = float(SeeStats.mean()(i, j));

            auto positionPair = GuiHelper::sphericalSurfacePositionPair(
                    electronsVector.positionsVector()[i], electronRadius,
                    electronsVector.positionsVector()[j], electronRadius);

            bool lineMode = false;
            if(lineMode) {
                if (std::abs(corr) >= spinCorrelationThreshold) {
                    if (corr < 0)
                        new Line3D(this, Qt::green, positionPair, std::abs(corr));
                    else
                        new Line3D(this, Qt::magenta, positionPair, std::abs(corr));
                }
            } else {
                if (std::abs(corr) >= spinCorrelationThreshold) {
                    Cylinder* c;
                    if (corr < 0) {
                        c = new Cylinder(this, Qt::green, positionPair, electronRadius / 7.5f, std::abs(corr));
                        c->material->setShininess(0);
                    }
                    //else
                    //    c = new Cylinder(this, Qt::magenta, positionPair, electronRadius / 7.5f, std::abs(corr));
                    //c->material->setShininess(0);
                }
            }

        }
    }
}
