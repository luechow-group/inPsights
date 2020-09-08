// Copyright (C) 2020 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

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
            const MolecularGeometry &permutee,
            const std::vector<Eigen::Index> &referenceIndicesInKitSystem,
            const MolecularGeometry &reference,
            double distanceMatrixCovarianceTolerance);
}


class EnvironmentBlock {
public:
    EnvironmentBlock(
            const std::deque<BestMatch::SOAPSimilarity::GrowingPerm>& possiblePerms,
            const MolecularGeometry &permutee,
            const MolecularGeometry &reference);

    void initialize(const std::deque<BestMatch::SOAPSimilarity::GrowingPerm>& possiblePerms);

    std::vector<std::vector<Eigen::Index>> filterPermutations(double distanceMatrixCovarianceTolerance);

    std::vector<std::vector<Eigen::Index>> permutedPermuteeIndicesCollection_; // of the permutee indices
    std::vector<Eigen::Index> permuteeIndices_, referenceIndices_;
    const MolecularGeometry &permutee_, &reference_;
};

class EnvironmentBlockJoiner {
public:
    EnvironmentBlockJoiner(const MolecularGeometry &permutee, const MolecularGeometry &reference);

    void initialize(const EnvironmentBlock &initialBlock);

    bool addBlock(const EnvironmentBlock &block, double distanceMatrixCovarianceTolerance);

    std::vector<std::vector<Eigen::Index>> jointPermutedPermuteeIndicesCollection_;
    std::vector<Eigen::Index> jointReferenceIndices_;
    const MolecularGeometry &permutee_, &reference_;
};


#endif //INPSIGHTS_ENVIRONMENTBLOCK_H
