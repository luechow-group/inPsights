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
                               + Elements::ElementInfo::vdwRadius(operator[](j).type()))/10.0;

            if ( distance < bondDrawingLimit)
                new Bond3D(*particles3D_[i], *particles3D_[j]);
        }
    }
}

template<>
void ElectronsVector3D::drawConnections() {

    double maxDistance = 1.6;
    double identicalThreshold = 0.01; //TODO add to general settings

    std::vector<long> pairedElectronsIndices, freeElectronsIndices;

    // find pairs
    for (long i = 0; i < numberOfEntities(); ++i) {
        bool paired = false;
        for (long j = i+1; j < numberOfEntities(); ++j) {
            if(Metrics::distance(operator[](i).position(),operator[](j).position()) <= identicalThreshold ) {
                paired = true;

                pairedElectronsIndices.push_back(i);
                pairedElectronsIndices.push_back(j);

                assert(operator[](i).type() != operator[](j).type()
                && "Antisymmetry violation! Types of paired electrons cannot be identical.");
            }

        }
        if (!paired)
            freeElectronsIndices.push_back(i);
    }

    std::sort(pairedElectronsIndices.begin(), pairedElectronsIndices.end());
    pairedElectronsIndices.erase(std::unique(
            pairedElectronsIndices.begin(),
            pairedElectronsIndices.end() ), pairedElectronsIndices.end()); // erase duplicates


    for (auto it = freeElectronsIndices.begin(); it != freeElectronsIndices.end(); ++it) {
            for (auto jt = it + 1; jt != freeElectronsIndices.end(); ++jt) {
                auto e1 = operator[](*it);
                auto e2 = operator[](*jt);

                auto q1 = GuiHelper::toQVector3D(e1.position());
                auto q2 = GuiHelper::toQVector3D(e2.position());

                if (Metrics::distance(e1.position(), e2.position()) < maxDistance) {
                    if (e1.type() == e2.type())
                        if (e1.type() == Spin::alpha)
                            Cylinder(connections_, GuiHelper::QColorFromType<Spin>(Spin::alpha), {q1, q2}, 0.015,
                                     0.5);//TODO connect to atoms
                        else if (e1.type() == Spin::beta)
                            Cylinder(connections_, GuiHelper::QColorFromType<Spin>(Spin::beta), {q1, q2}, 0.015,
                                     0.5);
                        else
                            Cylinder(connections_, GuiHelper::QColorFromType<Spin>(Spin::none), {q1, q2}, 0.015,
                                     0.5);
                    else
                        Line3D(connections_, Qt::black, {q1, q2}, 0.25);
                }
            }
        }

}