/* Copyright (C) 2020 Michael Heuer.
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

#ifndef INPSIGHTS_ENVIRONMENTBLOCK_H
#define INPSIGHTS_ENVIRONMENTBLOCK_H


#include <Combinatorics.h>
#include <Enumerate.h>
#include <MolecularGeometry.h>
#include <deque>
#include <BestMatchSimilarity.h>
#include <spdlog/spdlog.h>


namespace DistanceCovariance {
    bool conservingQ(
            const std::vector<Eigen::Index> &permuteeIndicesInKitSystem,
            const ElectronsVector &permutee,
            const std::vector<Eigen::Index> &referenceIndicesInKitSystem,
            const ElectronsVector &reference,
            double distanceMatrixCovarianceTolerance);
}


class EnvironmentBlock {
public:
    EnvironmentBlock(
            const std::deque<BestMatch::SOAPSimilarity::GrowingPerm>& possiblePerms,
            const ElectronsVector &permutee,
            const ElectronsVector &reference);

    void initialize(const std::deque<BestMatch::SOAPSimilarity::GrowingPerm>& possiblePerms);

    std::vector<std::vector<Eigen::Index>> filterPermutations(double distanceMatrixCovarianceTolerance);

    std::vector<std::vector<Eigen::Index>> permutedPermuteeIndicesCollection_; // of the permutee indices
    std::vector<Eigen::Index> permuteeIndices_, referenceIndices_;
    const ElectronsVector &permutee_, &reference_;
};

class EnvironmentBlockJoiner {
public:
    EnvironmentBlockJoiner(const ElectronsVector &permutee, const ElectronsVector &reference);

    void initialize(const EnvironmentBlock &initialBlock);

    bool addBlock(const EnvironmentBlock &block, double distanceMatrixCovarianceTolerance);

    std::vector<std::vector<Eigen::Index>> jointPermutedPermuteeIndicesCollection_;
    std::vector<Eigen::Index> jointReferenceIndices_;
    const ElectronsVector &permutee_, &reference_;
};


#endif //INPSIGHTS_ENVIRONMENTBLOCK_H
