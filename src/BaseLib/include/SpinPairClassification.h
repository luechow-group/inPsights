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

#ifndef INPSIGHTS_SPINPAIRCLASSIFICATION_H
#define INPSIGHTS_SPINPAIRCLASSIFICATION_H

#include <ParticlesVector.h>

namespace SpinPairClassification{
    enum class PairType { atSamePosition, closeBy, distant };

    using PairTypeMap = std::map<std::pair<Eigen::Index, Eigen::Index>, SpinPairClassification::PairType>;

    PairTypeMap classify(const ElectronsVector& electronsVector,
            double maxDistance = std::numeric_limits<double>::max(),
            double identicalThreshold = 0.01);

    bool atSamePositionQ(const std::pair<Eigen::Index, Eigen::Index>& pair,
            const ElectronsVector& electrons, double identicalThreshold = 0.01);

    bool closeByQ(const std::pair<Eigen::Index, Eigen::Index>& pair,
                         const ElectronsVector& electronsVector, double maxDistance = 1.6);

    bool isAtSamePositionQ(const PairTypeMap& pairTypes , Eigen::Index i);
}

#endif //INPSIGHTS_SPINPAIRCLASSIFICATION_H
