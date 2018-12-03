//
// Created by heuer on 03.12.18.
//

#include <SpinConnections3D.h>
#include <Metrics.h>
#include <Cylinder.h>
#include <Line3D.h>

enum class PairType {
    atSamePosition, closeBy
};

SpinConnections3D::SpinConnections3D(ElectronsVector3D *electronsVector3D, double maxDistance, double identicalThreshold)
        : IConnection(electronsVector3D->connections_),
        maxDistance_(maxDistance),
        identicalThreshold_(identicalThreshold),
        atSamePositionQList_(size_t(electronsVector3D->numberOfEntities()), false),
        pairIndicesMap_() {

    classifyElectrons(electronsVector3D);
    createConnections(electronsVector3D);
}

void SpinConnections3D::classifyElectrons(ElectronsVector3D *electronsVector3D) {
    for (long i = 0; i < electronsVector3D->numberOfEntities(); ++i) {
        for (long j = i + 1; j < electronsVector3D->numberOfEntities(); ++j) {
            auto dist = Metrics::distance(
                    electronsVector3D->operator[](i).position(),
                    electronsVector3D->operator[](j).position());

            if (dist <= identicalThreshold_) {
                atSamePositionQList_[i] = true;
                atSamePositionQList_[j] = true;
                pairIndicesMap_[{i, j}] = PairType::atSamePosition;

                assert(electronsVector3D->operator[](i).type() != electronsVector3D->operator[](j).type()
                       && "Antisymmetry violation! Types of paired electrons cannot be identical.");

            } else if (dist < maxDistance_) {
                pairIndicesMap_[{i, j}] = PairType::closeBy;
            }
        }
    }
}

void SpinConnections3D::createConnections(ElectronsVector3D *electronsVector3D) {
    for (auto &idxPair : pairIndicesMap_) {
        auto i = idxPair.first.first;
        auto j = idxPair.first.second;

        if (idxPair.second == PairType::closeBy && !atSamePositionQList_[i] && !atSamePositionQList_[j]) {

            auto e1 = electronsVector3D->operator[](i);
            auto e2 = electronsVector3D->operator[](j);

            auto q1 = GuiHelper::toQVector3D(e1.position());
            auto q2 = GuiHelper::toQVector3D(e2.position());

            if (Metrics::distance(e1.position(), e2.position()) < maxDistance_) {
                if (e1.type() == e2.type())
                    if (e1.type() == Spin::alpha)
                        new Cylinder(electronsVector3D->connections_,
                                     GuiHelper::QColorFromType<Spin>(Spin::alpha), {q1, q2}, 0.015, 0.5);
                    else if (e1.type() == Spin::beta)
                        new Cylinder(electronsVector3D->connections_,
                                     GuiHelper::QColorFromType<Spin>(Spin::beta), {q1, q2}, 0.015, 0.5);
                    else
                        new Cylinder(electronsVector3D->connections_,
                                     GuiHelper::QColorFromType<Spin>(Spin::none), {q1, q2}, 0.015, 0.5);
                else
                    new Line3D(electronsVector3D->connections_, Qt::black, {q1, q2}, 0.25);
            }
        }
    }
}
