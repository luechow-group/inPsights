//
// Created by Michael Heuer on 02.12.18.
//

#include <ParticlesVector3D.h>
#include <Metrics.h>
#include <Cylinder.h>
#include <Line3D.h>
#include <Bond3D.h>

template<>
void AtomsVector3D::drawConnections() {
    auto bondDrawingLimit = float(1.40 * 1e-10 / AU::length); //arbitrary choosen

    for (long i = 0; i < numberOfEntities(); ++i) {
        for (long j = 0; j < numberOfEntities(); ++j) {
            auto distance = Metrics::distance(operator[](i).position(), operator[](j).position())
                            - (Elements::ElementInfo::vdwRadius(operator[](i).type())
                               + Elements::ElementInfo::vdwRadius(operator[](j).type())) / 10.0;

            if (distance < bondDrawingLimit)
                new Bond3D(connections_, *particles3D_[i], *particles3D_[j]);
        }
    }
}

enum class PairType {
    atSamePosition, closeBy
};

template<>
void ElectronsVector3D::drawConnections() {

    double maxDistance = 1.6;
    double identicalThreshold = 0.01; //TODO add to general settings

    std::vector<bool> atSamePositionQList(numberOfEntities(), false);
    std::map<std::pair<long, long>, PairType> pairIndicesMap;

    // find pairs
    for (long i = 0; i < numberOfEntities(); ++i) {
        for (long j = i + 1; j < numberOfEntities(); ++j) {
            auto dist = Metrics::distance(operator[](i).position(), operator[](j).position());
            if (dist <= identicalThreshold) {

                atSamePositionQList[i] = true;
                atSamePositionQList[j] = true;
                pairIndicesMap[{i, j}] = PairType::atSamePosition;

                assert(operator[](i).type() != operator[](j).type()
                       && "Antisymmetry violation! Types of paired electrons cannot be identical.");

            } else if (dist < maxDistance) {
                pairIndicesMap[{i, j}] = PairType::closeBy;
            }
        }
    }

    for (auto &idxPair : pairIndicesMap) {
        auto i = idxPair.first.first;
        auto j = idxPair.first.second;

        if (idxPair.second == PairType::closeBy && !atSamePositionQList[i] && !atSamePositionQList[j]) {

            auto e1 = operator[](i);
            auto e2 = operator[](j);

            auto q1 = GuiHelper::toQVector3D(e1.position());
            auto q2 = GuiHelper::toQVector3D(e2.position());

            if (Metrics::distance(e1.position(), e2.position()) < maxDistance) {
                if (e1.type() == e2.type())
                    if (e1.type() == Spin::alpha)
                        new Cylinder(connections_, GuiHelper::QColorFromType<Spin>(Spin::alpha), {q1, q2}, 0.015, 0.5);
                    else if (e1.type() == Spin::beta)
                        new Cylinder(connections_, GuiHelper::QColorFromType<Spin>(Spin::beta), {q1, q2}, 0.015, 0.5);
                    else
                        new Cylinder(connections_, GuiHelper::QColorFromType<Spin>(Spin::none), {q1, q2}, 0.015, 0.5);
                else
                    new Line3D(connections_, Qt::black, {q1, q2}, 0.25);
            }
        }
    }

}