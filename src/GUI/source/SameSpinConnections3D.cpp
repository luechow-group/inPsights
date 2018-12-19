//
// Created by Michael Heuer on 04.12.18.
//

#include <SameSpinConnections3D.h>
#include <SpinClassification.h>
#include <Metrics.h>
#include <Cylinder.h>
#include <Line3D.h>

SameSpinConnections3D::SameSpinConnections3D(
        ElectronsVector3D *electronsVector3D,
        double maxDistance,
        double identicalThreshold)
        : IConnection(electronsVector3D->connections_) {

    createConnections(electronsVector3D, maxDistance, identicalThreshold);
};


void SameSpinConnections3D::createConnections(ElectronsVector3D *electronsVector3D, double maxDistance, double identicalThreshold) {

    auto pairTypes = SpinClassification::classify(*electronsVector3D, maxDistance, identicalThreshold );

    for (auto &idxPair : pairTypes) {
        auto i = idxPair.first.first;
        auto j = idxPair.first.second;

        if (idxPair.second == SpinClassification::PairType::closeBy
            && !SpinClassification::isAtSamePositionQ(pairTypes,i)
            && !SpinClassification::isAtSamePositionQ(pairTypes,j)
            ) {

            auto e1 = electronsVector3D->operator[](i);
            auto e2 = electronsVector3D->operator[](j);

            auto q1 = GuiHelper::toQVector3D(e1.position());
            auto q2 = GuiHelper::toQVector3D(e2.position());

            if (Metrics::distance(e1.position(), e2.position()) < maxDistance) {
                if (e1.type() == e2.type())
                    if (e1.type() == Spin::alpha)
                        new Cylinder(this,
                                     GuiHelper::QColorFromType<Spin>(Spin::alpha), {q1, q2}, 0.015, 0.5);
                    else if (e1.type() == Spin::beta)
                        new Cylinder(this,
                                     GuiHelper::QColorFromType<Spin>(Spin::beta), {q1, q2}, 0.015, 0.5);
                    else
                        new Cylinder(this,
                                     GuiHelper::QColorFromType<Spin>(Spin::none), {q1, q2}, 0.015, 0.5);
                else
                    new Line3D(this, Qt::black, {q1, q2}, 0.25);
            }
        }
    }
}