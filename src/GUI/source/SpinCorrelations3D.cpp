//
// Created by heuer on 03.12.18.
//

#include <SpinCorrelations3D.h>
#include <Line3D.h>

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

    auto sphereRadius = GuiHelper::radiusFromType(Spin::alpha);

    for (auto &idxPair : pairTypes) {
        auto i = idxPair.first.first;
        auto j = idxPair.first.second;

        if (idxPair.second == SpinPairClassification::PairType::closeBy
        && !SpinPairClassification::isAtSamePositionQ(pairTypes,i)
        && !SpinPairClassification::isAtSamePositionQ(pairTypes,j)
        ) {

            auto corr = SeeStats.mean()(i, j);
            if (std::abs(corr) >= spinCorrelationThreshold) {

                QColor color;
                if (corr > 0)
                    color = QColor::fromRgb(255, 0, 255);
                else
                    color = QColor::fromRgb(0, 255, 0);


                new Line3D(this, color, GuiHelper::sphericalSurfacePositionPair(
                        electronsVector.positionsVector()[i], sphereRadius,
                        electronsVector.positionsVector()[j], sphereRadius),
                                std::abs(corr));
            }
        }
    }
}
