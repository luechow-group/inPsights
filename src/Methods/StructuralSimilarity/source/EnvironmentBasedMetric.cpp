// Copyright (C) 2019-2020 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include "EnvironmentBasedMetric.h"
#include <StructuralSimilarity.h>
#include <Hungarian.h>
#include <limits>
#include <deque>
#include <utility>
#include <vector>
#include <Eigen/Core>
#include <ToString.h>
#include <EnvironmentBlock.h>
#include "GraphAnalysis.h"

/*
Eigen::MatrixXd Metrics::Similarity::EnvironmentBased::calculateEnvironmentSimilarityMatrix(
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
}*/

/*
 * Returns a list of distance-preserving permutations matching equivalent environments within a given
 * similarity threshold. If no match is found above the given threshold (within a comparison epsilon compensating
 * numerical imprecision), the best-match permutation and the most deviating environmental similarity along the
 * diagonal of the best-match permuted environmental similarity matrix is returned.
 *
 * The algorithm has five steps:
 *  1. Find matching dependent environments.
 *  2. Find possible permutations between dependent environments and check if they are distance preserving.
 *  3. Check for all permutations within these blocks for distance conservation.
 *  4. Combine the blockwise distance-conserving permutations to find overall distance-conserving permutations.
 *  5. Lastly, the distance-conserving permutations are converted back into the lab system
 *
 * Distance conservation is checked via the distance-covariance matrix difference between the permuted and unpermuted
 * permutee and the permuted and unpermuted reference.
 */
std::vector<BestMatch::DescendingMetricResult> Metrics::Similarity::EnvironmentBased::getBestMatchResults(
        const ::SOAP::MolecularSpectrum &permutee,
        const ::SOAP::MolecularSpectrum &reference,
        double distanceMatrixCovarianceTolerance,
        double similarityThreshold,
        double comparisionEpsilon) {

    assert(::SOAP::ParticleKit::isSubsetQ(permutee.molecule_)
           && "The permutee must be a subset of the particle kit.");
    assert(::SOAP::ParticleKit::isSubsetQ(reference.molecule_)
           && "The reference must be a subset of the particle kit.");

    const auto permuteeToKit = ::SOAP::ParticleKit::toKitPermutation(permutee.molecule_);
    const auto referenceToKit = ::SOAP::ParticleKit::fromKitPermutation(reference.molecule_);
    const auto referenceFromKit = ::SOAP::ParticleKit::fromKitPermutation(reference.molecule_);
    Eigen::PermutationMatrix<Eigen::Dynamic> identity(permutee.molecule_.numberOfEntities());

    spdlog::debug("Electrons of permutee {} in particle-kit system.",
            permuteeToKit.indices() == identity.indices() ?  "are": "are NOT");
    spdlog::debug("Electrons of reference {} in particle-kit system.",
                  referenceToKit.indices() == identity.indices() ?  "are": "are NOT");

    const auto particleNumber = size_t(permutee.molecule_.numberOfEntities());

    assert(permutee.molecule_.electrons().typesVector().countOccurence(Spin::alpha)
           == reference.molecule_.electrons().typesVector().countOccurence(Spin::alpha)
           && "The number of alpha electrons must match.");
    assert(permutee.molecule_.electrons().typesVector().countOccurence(Spin::beta)
           == reference.molecule_.electrons().typesVector().countOccurence(Spin::beta)
           && "The number of beta electrons must match.");
    assert(permutee.molecule_.atoms() == reference.molecule_.atoms()
           && "The atom vectors must match.");

    // environment similarity matrix of the permutee (rows) and the reference (cols)
    Eigen::MatrixXd environmentSimilarities = SOAP::StructuralSimilarity::correlationMatrix(permutee, reference);

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
    auto matches = Metrics::Similarity::EnvironmentBased::findEnvironmentMatches(environmentSimilarities, similarityThreshold, comparisionEpsilon);
    auto dependentMatches = Metrics::Similarity::EnvironmentBased::clusterDependentMatches(matches);
    for(const auto& dependentMatch : dependentMatches){
        if(dependentMatch.size() > 12)
            spdlog::warn("More than 12 equivalent environments found. This might indicate too indistinct SOAP settings.");
    }

    // Step 2.
    std::vector<EnvironmentBlock> blocks;
    for (const auto &dependentMatch : dependentMatches) {

        auto possiblePermutations = Metrics::Similarity::EnvironmentBased::findPossiblePermutations(dependentMatch);

        auto block = EnvironmentBlock(possiblePermutations, permutee.molecule_, reference.molecule_);

        auto filteredPerms = block.filterPermutations(distanceMatrixCovarianceTolerance);

        if (filteredPerms.empty()) { // no distant preserving permutation could be found
            spdlog::debug("exited early (not intra-block distance preserving");
            return {earlyExitResult};
        } else {
            blocks.emplace_back(block);
        }
    }

    // Step 3.
    EnvironmentBlockJoiner jointBlocks(permutee.molecule_, reference.molecule_);
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
        assert(jointPermutedPermuteeIndices.size() == size_t(particleNumber) &&
               "The found index ordering size must match the number of particles.");

        auto jointReferenceIndices = jointBlocks.jointReferenceIndices_;

        // construct permutations in the kit system
        Eigen::VectorXi finalPermIndicesInKitSystem(particleNumber);

        // determine final permutation indices and lowest environment similarity value for each electron
        double lowestEnvironmentSimilarityOfParticle = 1.0; // start with maximal value
        for (size_t i = 0; i < particleNumber; ++i) {

            auto permIdx = jointPermutedPermuteeIndices[i];
            auto refIdx = jointReferenceIndices[i];

            // create permutation vector
            finalPermIndicesInKitSystem[permIdx] = refIdx;
            spdlog::debug("({} {})", permIdx, refIdx);

            // determine lowest similarity value
            auto environmentSimilarity = environmentSimilarities(permIdx, refIdx);
            if (environmentSimilarity < lowestEnvironmentSimilarityOfParticle)
                lowestEnvironmentSimilarityOfParticle = environmentSimilarity;
        }
        spdlog::debug("Permutation (particle-kit system): {}", ToString::vectorXiToString(finalPermIndicesInKitSystem));

        Eigen::PermutationMatrix<Eigen::Dynamic> perm(finalPermIndicesInKitSystem);

        // Step 5.
        results.emplace_back(BestMatch::DescendingMetricResult({
            lowestEnvironmentSimilarityOfParticle,
            referenceFromKit * perm * permuteeToKit
        }));
    }
    std::sort(results.begin(), results.end()); // higher metric values come first

    for (auto r : results)
        spdlog::debug("Metric: {:01.16f}, Permutation: {} (lab system)", r.metric,
                      ToString::vectorXiToString(r.permutation.indices()));
    return results;
}


double Metrics::Similarity::EnvironmentBased::earlyExitMetric(const Eigen::MatrixXd &bestMatchPermutedEnvironmentSimilarities) {
    return bestMatchPermutedEnvironmentSimilarities.diagonal().minCoeff();
}

BestMatch::DescendingMetricResult Metrics::Similarity::EnvironmentBased::compare(
        const SOAP::MolecularSpectrum &permutee,
        const SOAP::MolecularSpectrum &reference,
        double distanceMatrixCovarianceTolerance,
        double soapThreshold,
        double numericalPrecisionEpsilon) {
    return Metrics::Similarity::EnvironmentBased::getBestMatchResults(permutee, reference, distanceMatrixCovarianceTolerance,
                                                                      soapThreshold, numericalPrecisionEpsilon).front();
}

Metrics::Similarity::EnvironmentBased::GrowingPerm::GrowingPerm(std::set<Eigen::Index> remainingPermuteeIndices,
                                                                std::deque<std::pair<Eigen::Index, Eigen::Index>> chainOfSwaps)
                                                    : remainingPermuteeIndices_(std::move(remainingPermuteeIndices)),
                                                      chainOfSwaps_(std::move(chainOfSwaps)){}


bool Metrics::Similarity::EnvironmentBased::GrowingPerm::add(const std::pair<Eigen::Index, Eigen::Index>& envMatch){

    auto permuteeIndexIterator = remainingPermuteeIndices_.find(envMatch.first);
    if(permuteeIndexIterator!= std::end(remainingPermuteeIndices_)) {
        remainingPermuteeIndices_.erase(permuteeIndexIterator);
        chainOfSwaps_.emplace_back(envMatch);
        return true;
    } else {
        return false;
    }
}

// Returns an object for each environment in the reference listing all equivalent environments in the permutee.
std::deque<Metrics::Similarity::EnvironmentBased::PermuteeToReferenceMatch>
Metrics::Similarity::EnvironmentBased::findEnvironmentMatches(
        const Eigen::MatrixXd &environmentSimilarities,
        double soapThreshold, double numericalPrecisionEpsilon) {

    auto adjacencyMatrix = GraphAnalysis::filter(environmentSimilarities, soapThreshold-numericalPrecisionEpsilon);

    std::deque<PermuteeToReferenceMatch> matches;
    for (Eigen::Index j = 0; j < environmentSimilarities.cols(); ++j) {
        auto equivalentPermuteeEnvs = GraphAnalysis::findVerticesOfIncomingEdges(adjacencyMatrix, j);
        matches.emplace_back(PermuteeToReferenceMatch{equivalentPermuteeEnvs, j});
    }

    return matches;
}

// finds PermuteeToReferenceMatch objects having overlapping permutee indices and groups them into lists.
std::deque<std::deque<Metrics::Similarity::EnvironmentBased::PermuteeToReferenceMatch>>
Metrics::Similarity::EnvironmentBased::clusterDependentMatches(const std::deque<PermuteeToReferenceMatch> &matches) {

    std::deque<std::deque<Metrics::Similarity::EnvironmentBased::PermuteeToReferenceMatch>> dependentMatchesGroups;
    assert(!matches.empty());

    for(const auto& match : matches){

        bool overlapBetweenPermuteeIndicesFoundQ = false;
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
                    overlapBetweenPermuteeIndicesFoundQ = true;
                    goto breakTwoLoops;
                }
            }
        }
        breakTwoLoops:
        if(!overlapBetweenPermuteeIndicesFoundQ)
            dependentMatchesGroups.emplace_back(std::deque<PermuteeToReferenceMatch>({match}));
    }

    return dependentMatchesGroups;
}

bool Metrics::Similarity::EnvironmentBased::GrowingPerm::operator<(const Metrics::Similarity::EnvironmentBased::GrowingPerm &rhs) const {
    assert(this->chainOfSwaps_.size() == rhs.chainOfSwaps_.size());
    for (size_t i = 0; i < rhs.chainOfSwaps_.size(); ++i)
        if(this->chainOfSwaps_[i] != rhs.chainOfSwaps_[i])
            return this->chainOfSwaps_[i] < rhs.chainOfSwaps_[i];

    return this->chainOfSwaps_[0] < rhs.chainOfSwaps_[0];
}

std::deque<Metrics::Similarity::EnvironmentBased::GrowingPerm> Metrics::Similarity::EnvironmentBased::findPossiblePermutations(
        const std::deque<Metrics::Similarity::EnvironmentBased::PermuteeToReferenceMatch>& dependentMatches){

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
Metrics::Similarity::EnvironmentBased::permuteIndicesFromKitSystem(const std::vector<Eigen::Index> &kitSystemIndices,
                                                                   const Eigen::PermutationMatrix<Eigen::Dynamic>& fromKitPermutation) {
    assert(size_t(fromKitPermutation.indices().size()) >= kitSystemIndices.size());

    std::vector<Eigen::Index> permutedIndices(kitSystemIndices.size());

    for (Eigen::size_t i = 0; i < kitSystemIndices.size(); ++i)
        permutedIndices[i] = fromKitPermutation.indices()[kitSystemIndices[i]];

    return permutedIndices;
}

/*
 * Calculates the positional distances of selected indices
 */
Eigen::MatrixXd
Metrics::Similarity::EnvironmentBased::calculateDistanceCovarianceMatrixOfSelectedIndices(const PositionsVector &positions,
                                                                                          const std::vector<Eigen::Index> &indices) {
    assert(size_t(positions.numberOfEntities()) >= indices.size());

    Eigen::MatrixXd positionalDistances(indices.size(), indices.size());

    for (Eigen::size_t i = 0; i < indices.size(); ++i)
        for (Eigen::size_t j = 0; j < indices.size(); ++j)
            positionalDistances(i, j) = (positions[indices[j]] -
                                         positions[indices[i]]).norm();

    return positionalDistances;
};