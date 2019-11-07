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

#include <SpinPairClassification.h>
#include <Metrics.h>

SpinPairClassification::PairTypeMap
SpinPairClassification::classify(const ElectronsVector& electronsVector, double maxDistance, double identicalThreshold) {

    PairTypeMap pairTypes;

    for (Eigen::Index i = 0; i < electronsVector.numberOfEntities(); ++i) {
        for (Eigen::Index j = i + 1; j < electronsVector.numberOfEntities(); ++j) {
            if (atSamePositionQ({i,j},electronsVector, identicalThreshold))
                pairTypes[{i, j}] = PairType::atSamePosition;

            else if (closeByQ({i,j}, electronsVector, maxDistance))
                pairTypes[{i, j}] = PairType::closeBy;
            else
                pairTypes[{i, j}] = PairType::distant;
        }
    }
    return pairTypes;
}

bool SpinPairClassification::atSamePositionQ(const std::pair<Eigen::Index, Eigen::Index>& pair,
                                         const ElectronsVector &electrons, double identicalThreshold) {
    auto dist = Metrics::distance(
            electrons[pair.first].position(),
            electrons[pair.second].position());

    if(dist <= identicalThreshold){
        assert(electrons[pair.first].type() != electrons[pair.second].type()
               && "Antisymmetry violation! Two same spin electrons cannot have the same position.");
        return true;
    } else {
        return false;
    }
}

bool SpinPairClassification::closeByQ(const std::pair<Eigen::Index, Eigen::Index>& pair,
              const ElectronsVector& electronsVector, double maxDistance) {
    auto dist = Metrics::distance(
            electronsVector[pair.first].position(),
            electronsVector[pair.second].position());

    return  dist <= maxDistance;
}


bool SpinPairClassification::isAtSamePositionQ(const PairTypeMap &pairTypes, Eigen::Index i) {

    for(const auto& pair : pairTypes)
        if( (pair.first.first == i || pair.first.second == i) && pair.second == PairType::atSamePosition)
            return true;
    return false;
}
