//
// Created by heuer on 03.12.18.
//

#include <SpinCorrelations3D.h>
#include <Line3D.h>

SpinCorrelations3D::SpinCorrelations3D(ElectronsVector3D *electronsVector3D,
                                       Statistics::RunningStatistics<Eigen::MatrixXd, unsigned, true> SeeStats,
                                       double spinCorrelationThreshold)
        :
        SpinConnections3D(electronsVector3D, 0.01),
        SeeStats_(std::move(SeeStats)) {

    drawSpinCorrelations(electronsVector3D, spinCorrelationThreshold);
}

void SpinCorrelations3D::drawSpinCorrelations(ElectronsVector3D *electronsVector3D, double spinCorrelationThreshold)  {

    for (auto &idxPair : pairIndicesMap_) {
        auto i = idxPair.first.first;
        auto j = idxPair.first.second;

        if (idxPair.second == PairType::closeBy && !atSamePositionQList_[i] && !atSamePositionQList_[j]) {

            auto corr = SeeStats_.mean()(i, j);
            if (std::abs(corr) >= spinCorrelationThreshold) {

                QColor color;
                if (corr > 0)
                    color = QColor::fromRgb(255, 0, 255);
                else
                    color = QColor::fromRgb(0, 255, 0);

                new Line3D(electronsVector3D->correlations_, color, {
                        GuiHelper::toQVector3D(electronsVector3D->positionsVector()[i]),
                        GuiHelper::toQVector3D(electronsVector3D->positionsVector()[j])}, std::abs(corr));
            }
        }
    }
}
