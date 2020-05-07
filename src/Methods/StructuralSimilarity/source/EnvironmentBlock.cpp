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

#include <EnvironmentBlock.h>

/*
 * Distance conservation is checked by calculating the difference of the distance covariance matrice of
 * the permutee and the reference, and assuring, that the maximal distance is below a certain threshold.
 */
bool DistanceCovariance::conservingQ(
        const std::vector<Eigen::Index> &permuteeIndices,
        const ElectronsVector &permutee,
        const std::vector<Eigen::Index> &referenceIndices,
        const ElectronsVector &reference,
        double distanceMatrixCovarianceTolerance) {

    auto covA = BestMatch::SOAPSimilarity::calculateDistanceCovarianceMatrixOfSelectedIndices(permutee, permuteeIndices);
    auto covB = BestMatch::SOAPSimilarity::calculateDistanceCovarianceMatrixOfSelectedIndices(reference, referenceIndices);

    auto conservingQ = (covB - covA).array().abs().maxCoeff() <= distanceMatrixCovarianceTolerance;

    if (!conservingQ)
        spdlog::debug("Distance covariance matrix difference:\n{}",
                ToString::matrixXdToString((covB - covA), 3));

    return conservingQ;
};

EnvironmentBlock::EnvironmentBlock(
        std::deque<std::pair<Eigen::Index, Eigen::Index>> indexPairs,
        const ElectronsVector &permutee,
        const ElectronsVector &reference)
        : permuteeIndices_(),
          referenceIndices_(),
          permutee_(permutee),
          reference_(reference) {
    initialize(indexPairs);
}


void EnvironmentBlock::initialize(std::deque<std::pair<Eigen::Index, Eigen::Index>> indexPairs) {
    for (const auto &indexPair : indexPairs) {
        permuteeIndices_.emplace_back(indexPair.first);
        referenceIndices_.emplace_back(indexPair.second);
    }
    std::sort(permuteeIndices_.begin(), permuteeIndices_.end());
    std::sort(referenceIndices_.begin(), referenceIndices_.end());

    permutedPermuteeIndicesCollection_ = Combinatorics::Permutations<Eigen::Index>(permuteeIndices_).get();

    spdlog::debug("Reference indices of current block:");
    spdlog::debug(ToString::stdvectorLongIntToString(referenceIndices_));
    spdlog::debug("Permuted permutee indices for current block:");

    for(const auto& permutedPermuteeIndices : permutedPermuteeIndicesCollection_)
        spdlog::debug(ToString::stdvectorLongIntToString(permutedPermuteeIndices));
};

std::vector<std::vector<Eigen::Index>> EnvironmentBlock::filterPermutations(double distanceMatrixCovarianceTolerance) {

    for (auto it = permutedPermuteeIndicesCollection_.begin(); it != permutedPermuteeIndicesCollection_.end(); it++) {
        auto conservingQ = DistanceCovariance::conservingQ(
                *it, permutee_,
                referenceIndices_, reference_,
                distanceMatrixCovarianceTolerance);
        if (!conservingQ)
            permutedPermuteeIndicesCollection_.erase(it--); // decrements iterator after erase
    }
    return permutedPermuteeIndicesCollection_;
};


EnvironmentBlockJoiner::EnvironmentBlockJoiner(const ElectronsVector &permutee, const ElectronsVector &reference)
        : jointPermutedPermuteeIndicesCollection_(0),
          jointReferenceIndices_(0),
          permutee_(permutee),
          reference_(reference) {}

void EnvironmentBlockJoiner::initialize(const EnvironmentBlock &initialBlock) {
    jointPermutedPermuteeIndicesCollection_ = initialBlock.permutedPermuteeIndicesCollection_;
    jointReferenceIndices_ = initialBlock.referenceIndices_;
}

bool EnvironmentBlockJoiner::addBlock(const EnvironmentBlock &block, double distanceMatrixCovarianceTolerance) {
    assert(block.permutee_ == permutee_);
    assert(block.reference_ == reference_);

    // initialize
    if (jointReferenceIndices_.empty() && jointPermutedPermuteeIndicesCollection_.empty()) {
        initialize(block);
        return true; // blocks are already distance conserving
    }

    // add reference indices of new block and don't permute them
    for (auto i : block.referenceIndices_) {
        jointReferenceIndices_.emplace_back(i);
    }

    spdlog::debug("Reference indices of current block:");
    spdlog::debug(ToString::stdvectorLongIntToString(jointReferenceIndices_));


    std::vector<std::vector<Eigen::Index>> updatedJointPermutedPermuteeIndicesCollection;

    // for each joint permutation that survived so far, append and test all permutations from the new block
    for (const auto &jointPermutedPermuteeIndices : jointPermutedPermuteeIndicesCollection_) {

        // insert new indices
        for (auto newIndices: block.permutedPermuteeIndicesCollection_) {
            auto jointPermutedPermuteeIndicesCopy = jointPermutedPermuteeIndices;

            // start with copy and add new indices
            jointPermutedPermuteeIndicesCopy.insert(
                    jointPermutedPermuteeIndicesCopy.end(),
                    newIndices.begin(), newIndices.end());

            spdlog::debug("Try new joint permuted permutee indices:");
            spdlog::debug(ToString::stdvectorLongIntToString(jointPermutedPermuteeIndicesCopy));

            assert(jointPermutedPermuteeIndicesCopy.size() == jointReferenceIndices_.size());
            auto conservingQ = DistanceCovariance::conservingQ(
                    jointPermutedPermuteeIndicesCopy, permutee_,
                    jointReferenceIndices_, reference_,
                    distanceMatrixCovarianceTolerance);

            spdlog::debug("{}",conservingQ? "OK" : "NOT CONSERVING");

            if (conservingQ)
                updatedJointPermutedPermuteeIndicesCollection.emplace_back(jointPermutedPermuteeIndicesCopy);
        }
    }

    jointPermutedPermuteeIndicesCollection_ = updatedJointPermutedPermuteeIndicesCollection;

    return !jointPermutedPermuteeIndicesCollection_.empty();
}