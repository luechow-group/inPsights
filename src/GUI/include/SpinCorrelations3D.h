//
// Created by heuer on 03.12.18.
//

#ifndef INPSIGHTS_SPINCORRELATIONS3D_H
#define INPSIGHTS_SPINCORRELATIONS3D_H

#include <IConnection.h>
#include <ParticlesVector3D.h>

#include <Statistics.h>
#include <Line3D.h>

class SpinCorrelations3D : public IConnection {
public:
    SpinCorrelations3D(ElectronsVector3D *electronsVector3D,
            const Statistics::RunningStatistics<Eigen::MatrixXd, unsigned, true> &SeeStats,
            int spinCorrelationThreshold)
            : IConnection(electronsVector3D->correlations_) {

        drawSpinCorrelations(electronsVector3D, SeeStats, spinCorrelationThreshold);
    }

    void drawSpinCorrelations(ElectronsVector3D *electronsVector3D,
                              const Statistics::RunningStatistics<Eigen::MatrixXd, unsigned, true> &SeeStats,
                              int spinCorrelationThreshold) {

        for (long i = 0; i < electronsVector3D->numberOfEntities(); ++i) {
            for (long j = i + 1; j < electronsVector3D->numberOfEntities(); ++j) {

                auto corr = SeeStats.mean()(i, j);
                if (int(std::abs(corr)*255) >= spinCorrelationThreshold) {

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
};

#endif //INPSIGHTS_SPINCORRELATIONS3D_H
