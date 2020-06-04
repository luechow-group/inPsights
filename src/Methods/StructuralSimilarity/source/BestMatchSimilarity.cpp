/* Copyright (C) 2019-2020 Michael Heuer.
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

#include "BestMatchSimilarity.h"
#include <LocalSimilarity.h>
#include <Hungarian.h>
#include <limits>
#include <deque>
#include <utility>
#include <vector>
#include <Eigen/Core>
#include <ToString.h>
#include <EnvironmentBlock.h>
#include <GraphAnalysis.h>

using namespace SOAP;

Eigen::MatrixXd BestMatch::SOAPSimilarity::calculateEnvironmentSimilarityMatrix(
        const MolecularSpectrum &permutee,
        const MolecularSpectrum &reference) {

    auto nAlpha = reference.molecule_.electrons().typesVector().countOccurence(Spin::alpha);
    auto nBeta = reference.molecule_.electrons().typesVector().countOccurence(Spin::beta);

    auto N = nAlpha + nBeta;
    Eigen::MatrixXd environmentalSimilarities(N, N);

    auto zeta = General::settings.zeta();

    // TODO consider identical spin flip for chemical mode?

    TypeSpecificNeighborhoodsAtOneCenter expA, expB;
    for (unsigned i = 0; i < nAlpha; ++i) {
        EnumeratedType<int> enumeratedType_i(Spins::spinToInt(Spin::alpha), i);
        expA = permutee.molecularCenters_.find(enumeratedType_i)->second;

        for (unsigned j = 0; j < nAlpha; ++j) {
            EnumeratedType<int> enumeratedType_j(Spins::spinToInt(Spin::alpha), j);
            expB = reference.molecularCenters_.find(enumeratedType_j)->second;
            environmentalSimilarities(i, j) = LocalSimilarity::kernel(expA, expB, zeta);
        }
        for (unsigned j = 0; j < nBeta; ++j) {
            EnumeratedType<int> enumeratedType_j(Spins::spinToInt(Spin::beta), j);
            expB = reference.molecularCenters_.find(enumeratedType_j)->second;
            environmentalSimilarities(i, nAlpha + j) = LocalSimilarity::kernel(expA, expB, zeta);
        }
    }
    for (unsigned i = 0; i < nBeta; ++i) {
        EnumeratedType<int> enumeratedType_i(Spins::spinToInt(Spin::beta), i);
        expA = permutee.molecularCenters_.find(enumeratedType_i)->second;

        for (unsigned j = 0; j < nAlpha; ++j) {
            EnumeratedType<int> enumeratedType_j(Spins::spinToInt(Spin::alpha), j);
            expB = reference.molecularCenters_.find(enumeratedType_j)->second;
            environmentalSimilarities(nAlpha + i, j) = LocalSimilarity::kernel(expA, expB, zeta);
        }
        for (unsigned j = 0; j < nBeta; ++j) {
            EnumeratedType<int> enumeratedType_j(Spins::spinToInt(Spin::beta), j);
            expB = reference.molecularCenters_.find(enumeratedType_j)->second;
            environmentalSimilarities(nAlpha + i, nAlpha + j) = LocalSimilarity::kernel(expA, expB, zeta);
        }
    }
    return environmentalSimilarities;
}

/*
 * Returns a list of distance-preserving permutations matching equivalent environments within a given
 * similarity threshold. If no match is found above the given threshold (within a comparision epsilion compensating
 * numerical imprecisions), the best-match permutation and the most deviating environmental similarity along the
 * diagonal of the best-match permuted environmental similarity matrix is returned.
 *
 * The algorithm has five steps:
 *  1. Find matching dependent environments.
 *  2. Find possible permutations between dependent environments and check of there are distance preserving.
 *  3. Check for all permutations within these blocks for distance conservation.
 *  4. Combine the blockwise distance-conserving permutations to find overall distance-conserving permutations.
 *  5. Lastly, the distance-conserving permutations are converted back into the lab system
 *
 * Distance conservation is checked via the distance-covariance matrix difference between the permuted and unpermuted
 * permutee and the permuted and unpermuted reference.
 */
std::vector<BestMatch::DescendingMetricResult> BestMatch::SOAPSimilarity::getBestMatchResults(
        const MolecularSpectrum &permutee,
        const MolecularSpectrum &reference,
        double distanceMatrixCovarianceTolerance,
        double similarityThreshold,
        double comparisionEpsilon) {

    assert(ParticleKit::isSubsetQ(permutee.molecule_)
           && "The permutee must be a subset of the particle kit.");
    assert(ParticleKit::isSubsetQ(reference.molecule_)
           && "The reference must be a subset of the particle kit.");

    const auto permuteeToKit = ParticleKit::toKitPermutation(permutee.molecule_.electrons());
    const auto referenceToKit = ParticleKit::fromKitPermutation(reference.molecule_.electrons());
    const auto referenceFromKit = ParticleKit::fromKitPermutation(reference.molecule_.electrons());
    Eigen::PermutationMatrix<Eigen::Dynamic> identity(permutee.molecule_.electrons().numberOfEntities());

    spdlog::debug("Electrons of permutee {} in particle-kit system.",
            permuteeToKit.indices() == identity.indices() ?  "are": "are NOT");
    spdlog::debug("Electrons of reference {} in particle-kit system.",
                  referenceToKit.indices() == identity.indices() ?  "are": "are NOT");

    const auto N = size_t(permutee.molecule_.electrons().numberOfEntities());

    assert(permutee.molecule_.electrons().typesVector().countOccurence(Spin::alpha)
           == reference.molecule_.electrons().typesVector().countOccurence(Spin::alpha)
           && "The number of alpha electrons must match.");
    assert(permutee.molecule_.electrons().typesVector().countOccurence(Spin::beta)
           == reference.molecule_.electrons().typesVector().countOccurence(Spin::beta)
           && "The number of beta electrons must match.");

    // environment similarity matrix of the permutee (rows) and the reference (cols)
    Eigen::MatrixXd environmentSimilarities = calculateEnvironmentSimilarityMatrix(permutee, reference);
    spdlog::debug("Environmental similarity matrix "
                  "(in the particle kit system, permutee = rows, reference = cols)\n{}",
                  ToString::matrixXdToString(environmentSimilarities, 3));


    auto bestMatch = Hungarian<double>::findMatching(environmentSimilarities, Matchtype::MAX);
    auto earlyExitResult = BestMatch::DescendingMetricResult({
        // find the best-match (row) permutation of the environments (of the permutee) maximizing the diagonal)
            earlyExitMetric(bestMatch * environmentSimilarities),
            referenceFromKit * bestMatch * permuteeToKit //convert best-match perm from particle-kit into the lab system
    });

    if (earlyExitResult.metric < (similarityThreshold - comparisionEpsilon)) {
        spdlog::debug("exited early (similarity < (threshold - epsilon) {:01.16f} < {:01.16f} - {:01.16f})",
                earlyExitResult.metric, similarityThreshold, comparisionEpsilon);

        return {earlyExitResult};
    }

    // Step 1.
    auto matches = BestMatch::SOAPSimilarity::findEnvironmentMatches(environmentSimilarities, similarityThreshold, comparisionEpsilon);
    auto dependentMatches = BestMatch::SOAPSimilarity::groupDependentMatches(matches);
    for(const auto& dependentMatch : dependentMatches){
        if(dependentMatch.size() > 12)
            spdlog::warn("More than 12 equivalent electronic environments found. This might indicate too indistinct SOAP settings.");
    }

    // Step 2.
    std::vector<EnvironmentBlock> blocks;
    for (const auto &dependentMatch : dependentMatches) {

        auto possiblePermutations = BestMatch::SOAPSimilarity::findPossiblePermutations(dependentMatch);

        auto block = EnvironmentBlock(possiblePermutations,
                                      permutee.molecule_.electrons(),
                                      reference.molecule_.electrons());

        auto filteredPerms = block.filterPermutations(distanceMatrixCovarianceTolerance);

        if (filteredPerms.empty()) { // no distant preserving permutation could be found
            spdlog::debug("exited early (not intra-block distance preserving");
            return {earlyExitResult};
        } else {
            blocks.emplace_back(block);
        }
    }

    // Step 3.
    EnvironmentBlockJoiner jointBlocks(permutee.molecule_.electrons(), reference.molecule_.electrons());
    for (const auto &block : blocks) {
        auto conservingQ = jointBlocks.addBlock(block, distanceMatrixCovarianceTolerance);

        if (!conservingQ) {
            spdlog::debug("exited early (not inter-block distance preserving)");
            return {earlyExitResult};
        }
    }

    // Step 4.
    std::vector<BestMatch::DescendingMetricResult> results;
    for (const auto &jointPermutedPermuteeIndices : jointBlocks.jointPermutedPermuteeIndicesCollection_) {
        assert(jointPermutedPermuteeIndices.size() == size_t(N) &&
               "The found index ordering size must match the number of electrons.");

        auto jointReferenceIndices = jointBlocks.jointReferenceIndices_;

        // construct permutations in the kit system
        Eigen::VectorXi finalPermIndicesInKitSystem(N);

        // determine final permutation indices and lowest environment similarity value for each electron
        double lowestEnvironmentSimilarityOfParticle = 1.0; // start with maximal value
        for (size_t i = 0; i < N; ++i) {

            auto permIdx = jointPermutedPermuteeIndices[i];
            auto refIdx = jointReferenceIndices[i];

            // create permuation vector
            finalPermIndicesInKitSystem[permIdx] = refIdx;
            spdlog::debug("({} {})", permIdx, refIdx);

            // determine lowest similarity
            auto environmentSimilarity = environmentSimilarities(permIdx, refIdx);
            if (environmentSimilarity < lowestEnvironmentSimilarityOfParticle)
                lowestEnvironmentSimilarityOfParticle = environmentSimilarity;
        }
        spdlog::debug("Permutation (particle-kit system): {}", ToString::vectorXiToString(finalPermIndicesInKitSystem));

        Eigen::PermutationMatrix<Eigen::Dynamic> perm(finalPermIndicesInKitSystem);

        results.emplace_back(BestMatch::DescendingMetricResult(
                {lowestEnvironmentSimilarityOfParticle, referenceFromKit * perm * permuteeToKit}));
    }
    std::sort(results.begin(), results.end()); // higher metric values come first

    for (auto r : results)
        spdlog::debug("Metric: {:01.16f}, Permutation: {} (lab system)", r.metric,
                      ToString::vectorXiToString(r.permutation.indices()));
    return results;
}


double BestMatch::SOAPSimilarity::earlyExitMetric(const Eigen::MatrixXd &bestMatchPermutedEnvironmentSimilarities) {
    return bestMatchPermutedEnvironmentSimilarities.diagonal().minCoeff();
}

BestMatch::DescendingMetricResult BestMatch::SOAPSimilarity::compare(
        const MolecularSpectrum &permutee,
        const MolecularSpectrum &reference,
        double distanceMatrixCovarianceTolerance,
        double soapThreshold,
        double numericalPrecisionEpsilon) {
    return BestMatch::SOAPSimilarity::getBestMatchResults(permutee, reference, distanceMatrixCovarianceTolerance,
                                                          soapThreshold, numericalPrecisionEpsilon).front();
}

BestMatch::SOAPSimilarity::GrowingPerm::GrowingPerm(std::set<Eigen::Index> remainingPermuteeIndices,
                                                    std::deque<std::pair<Eigen::Index, Eigen::Index>> chainOfSwaps)
                                                    : remainingPermuteeIndices_(std::move(remainingPermuteeIndices)),
                                                      chainOfSwaps_(std::move(chainOfSwaps)){}


bool BestMatch::SOAPSimilarity::GrowingPerm::add(const std::pair<Eigen::Index, Eigen::Index>& envMatch){

    auto permuteeIndexIterator = remainingPermuteeIndices_.find(envMatch.first);
    if(permuteeIndexIterator!= std::end(remainingPermuteeIndices_)) {
        remainingPermuteeIndices_.erase(permuteeIndexIterator);
        chainOfSwaps_.emplace_back(envMatch);
        return true;
    } else {
        return false;
    }
}

/*
 * Returns list of environments in the permutee matching a reference.
 * Multiple environments can be found for each reference index.
 */
std::deque<BestMatch::SOAPSimilarity::PermuteeToReferenceMatch>
BestMatch::SOAPSimilarity::findEnvironmentMatches(
        const Eigen::MatrixXd &environmentSimilarities,
        double soapThreshold, double numericalPrecisionEpsilon) {

    auto adjacencyMatrix = GraphAnalysis::filter(environmentSimilarities, soapThreshold-numericalPrecisionEpsilon);

    std::deque<PermuteeToReferenceMatch> matches;
    for (Eigen::Index j = 0; j < environmentSimilarities.cols(); ++j) {
        auto permuteeEnvs = GraphAnalysis::findVerticesOfIncomingEdges(adjacencyMatrix, j);
        matches.emplace_back(PermuteeToReferenceMatch{permuteeEnvs, j});
    }

    return matches;
}

std::deque<std::deque<BestMatch::SOAPSimilarity::PermuteeToReferenceMatch>>
BestMatch::SOAPSimilarity::groupDependentMatches(const std::deque<PermuteeToReferenceMatch> &matches) {

    std::deque<std::deque<BestMatch::SOAPSimilarity::PermuteeToReferenceMatch>> dependentMatchesGroups;
    assert(!matches.empty());

    for(const auto& match : matches){

        bool foundQ = false;
        for ( auto& dependentMatchesGroup : dependentMatchesGroups){

            for ( auto dependentMatch : dependentMatchesGroup ) {
                // check if index intersection between a match and any of the dependent matches exists
                std::list<Eigen::Index> intersection;
                std::set_intersection(
                        std::begin(match.permuteeIndices), std::end(match.permuteeIndices),
                        std::begin(dependentMatch.permuteeIndices), std::end(dependentMatch.permuteeIndices),
                        std::inserter(intersection,std::begin(intersection)));

                if(!intersection.empty()) {
                    dependentMatchesGroup.emplace_back(match);
                    foundQ = true;
                    goto breakTwoLoops;
                }
            }
        }
        breakTwoLoops:
        if(!foundQ)
            dependentMatchesGroups.emplace_back(std::deque<PermuteeToReferenceMatch>({match}));
    }

    return dependentMatchesGroups;
}

bool BestMatch::SOAPSimilarity::GrowingPerm::operator<(const BestMatch::SOAPSimilarity::GrowingPerm &rhs) const {
    assert((*this).chainOfSwaps_.size() == rhs.chainOfSwaps_.size());
    for (size_t i = 0; i < rhs.chainOfSwaps_.size(); ++i)
        if((*this).chainOfSwaps_[i] != rhs.chainOfSwaps_[i])
            return (*this).chainOfSwaps_[i] < rhs.chainOfSwaps_[i];

    return (*this).chainOfSwaps_[0] < rhs.chainOfSwaps_[0];
}

std::deque<BestMatch::SOAPSimilarity::GrowingPerm> BestMatch::SOAPSimilarity::findPossiblePermutations(
        const std::deque<BestMatch::SOAPSimilarity::PermuteeToReferenceMatch>& dependentMatches){

    assert(!dependentMatches.empty());

    std::set<Eigen::Index> allIndices;
    for(auto i : dependentMatches)
        allIndices.insert(std::begin(i.permuteeIndices), std::end(i.permuteeIndices));

    spdlog::debug("Number of dependent indices: {}", allIndices.size());

    std::deque<GrowingPerm> possiblePerms = {{allIndices,{}}};

    // try new perms by adding all swaps for all dependent matches
    for(auto& dependentMatch : dependentMatches){
        spdlog::debug("Number of possible perms: {}", possiblePerms.size());
        std::deque<GrowingPerm> newPossiblePerms;
        for(auto permuteeEnvIndex : dependentMatch.permuteeIndices){

            for (auto perm : possiblePerms) {
                auto possibleQ = perm.add({permuteeEnvIndex, dependentMatch.referenceIndex});
                if (possibleQ)
                    newPossiblePerms.emplace_back(perm);
            }
        }
        if(newPossiblePerms.empty())
            return {};
        else
            possiblePerms = newPossiblePerms;
    }

    std::sort(possiblePerms.begin(), possiblePerms.end());

    return possiblePerms;
};

std::vector<Eigen::Index>
BestMatch::SOAPSimilarity::permuteIndicesFromKitSystem(const std::vector<Eigen::Index> &kitSystemIndices,
                                                       const Eigen::PermutationMatrix<Eigen::Dynamic>& fromKitPermutation) {
    assert(fromKitPermutation.indices().size() >= kitSystemIndices.size());

    std::vector<Eigen::Index> permutedIndices(kitSystemIndices.size());

    for (Eigen::size_t i = 0; i < kitSystemIndices.size(); ++i)
        permutedIndices[i] = fromKitPermutation.indices()[kitSystemIndices[i]];

    return permutedIndices;
}

/*
 * Calculates the positional distances of selected indices
 */
Eigen::MatrixXd
BestMatch::SOAPSimilarity::calculateDistanceCovarianceMatrixOfSelectedIndices(const PositionsVector &positions,
                                                                              const std::vector<Eigen::Index> &indices) {
    assert(size_t(positions.numberOfEntities()) >= indices.size());

    Eigen::MatrixXd positionalDistances(indices.size(), indices.size());

    for (Eigen::size_t i = 0; i < indices.size(); ++i)
        for (Eigen::size_t j = 0; j < indices.size(); ++j)
            positionalDistances(i, j) = (positions[indices[j]] -
                                         positions[indices[i]]).norm();

    return positionalDistances;
};