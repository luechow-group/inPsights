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

#include <SameSpinConnections3D.h>
#include <SpinPairClassification.h>
#include <Metrics.h>
#include <Cylinder.h>
#include <Line3D.h>

SameSpinConnections3D::SameSpinConnections3D(
        ElectronsVector3D *electronsVector3D,
        double maxDistance,
        double identicalThreshold)
        : IConnection(electronsVector3D->connections_) {

    createConnections(*electronsVector3D, maxDistance, identicalThreshold);
};


void SameSpinConnections3D::createConnections(const ElectronsVector &electronsVector, double maxDistance,
                                              double identicalThreshold) {

    auto pairTypes = SpinPairClassification::classify(electronsVector, maxDistance, identicalThreshold );

    for (auto &idxPair : pairTypes) {
        auto i = idxPair.first.first;
        auto j = idxPair.first.second;

        if (idxPair.second == SpinPairClassification::PairType::closeBy
            && !SpinPairClassification::isAtSamePositionQ(pairTypes,i)
            && !SpinPairClassification::isAtSamePositionQ(pairTypes,j)
            ) {

            auto e1 = electronsVector[i];
            auto e2 = electronsVector[j];

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