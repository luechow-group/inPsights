/* Copyright (C) 2019 Michael Heuer.
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
#include <set>
#include <vector>
#include <Eigen/Core>
#include <ToString.h>
#include <EnvironmentBlock.h>

using namespace SOAP;

Eigen::MatrixXd BestMatch::SOAPSimilarity::calculateEnvironmentalSimilarityMatrix(
        const MolecularSpectrum &permutee,
        const MolecularSpectrum &reference) {

    auto nAlpha = reference.molecule_.electrons().typesVector().countOccurence(Spin::alpha);
    auto nBeta = reference.molecule_.electrons().typesVector().countOccurence(Spin::beta);

    auto N = nAlpha + nBeta;
    Eigen::MatrixXd environmentalSimilarities(N, N);

    auto zeta = General::settings.zeta();

    // TODO consider identical spin flip?

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

std::vector<BestMatch::DescendingMetricResult> BestMatch::SOAPSimilarity::getBestMatchResults(
        const MolecularSpectrum &permutee,
        const MolecularSpectrum &reference,
        double distanceMatrixCovarianceTolerance, double soapThreshold) {

    assert(ParticleKit::isSubsetQ(permutee.molecule_)
           && "The permutee must be a subset of the particle kit.");
    assert(ParticleKit::isSubsetQ(reference.molecule_)
           && "The reference must be a subset of the particle kit.");

    auto permuteeToKit = ParticleKit::toKitPermutation(permutee.molecule_.electrons());
    auto permuteeFromKit = ParticleKit::toKitPermutation(permutee.molecule_.electrons());
    auto referenceToKit = ParticleKit::fromKitPermutation(reference.molecule_.electrons());
    auto referenceFromKit = ParticleKit::fromKitPermutation(reference.molecule_.electrons());
    Eigen::PermutationMatrix<Eigen::Dynamic> identity(permutee.molecule_.electrons().numberOfEntities());

    spdlog::debug("Electrons of permutee are{} in particle-kit system.",
            permuteeToKit.indices() == identity.indices() ?  "": " NOT");
    spdlog::debug("Electrons of reference are{} in particle-kit system.",
                  referenceToKit.indices() == identity.indices() ?  "": " NOT");


    // TODO assert that identical number of electrons and same atom geometry? Is this constraint needed? What happens with rows/cols of zero?
    auto N = size_t(permutee.molecule_.electrons().numberOfEntities());

    assert(permutee.molecule_.electrons().typesVector().countOccurence(Spin::alpha)
           == reference.molecule_.electrons().typesVector().countOccurence(Spin::alpha)
           && "The number of alpha electrons must match.");
    assert(permutee.molecule_.electrons().typesVector().countOccurence(Spin::beta)
           == reference.molecule_.electrons().typesVector().countOccurence(Spin::beta)
           && "The number of beta electrons must match.");

    // environment matrix of the permutee (rows) and the reference (cols)
    Eigen::MatrixXd environmentalSimilarities = calculateEnvironmentalSimilarityMatrix(permutee, reference);
    spdlog::debug("Environmental similarity matrix "
                  "(in the particle kit system, permutee = rows, reference = cols)\n{}",
                  ToString::matrixXdToString(environmentalSimilarities, 3));

    // find the best match permutation of the environments (of the permutee => rows are permuted) that maximizes the diagonal)
    auto bestMatch = Hungarian<double>::findMatching(environmentalSimilarities, Matchtype::MAX);

    auto bestMatchPermutedEnvironmentalSimilarities = bestMatch * environmentalSimilarities;

    spdlog::debug("Best-Match permutation: {}\n Best-match permuted environmental similarity matrix "
                  "(in the particle kit system, permutee = rows, reference = cols)\n{}",
                  ToString::vectorXiToString(bestMatch.indices()),
                  ToString::matrixXdToString(bestMatchPermutedEnvironmentalSimilarities, 3));

    auto earlyExitResult = BestMatch::DescendingMetricResult(
            {earlyExitMetric(bestMatchPermutedEnvironmentalSimilarities), bestMatch});

    if (earlyExitResult.metric < soapThreshold) {
        spdlog::debug("exited early (similarity below threshold)");
        return {earlyExitResult};
    }

    // some indices might depend on each other (are equivalent, e.g. two electrons in a nucleus might be swapped)
    // obtain dependent indices in the kit system (blocks of indices might )
    auto blockwiseDependentIndexPairs =
            BestMatch::SOAPSimilarity::getBlockwiseDependentIndexPairs(
                    environmentalSimilarities,
                    bestMatch,soapThreshold);

    spdlog::debug("Dependent indice pairs");
    for (const auto &dependentIndexPairs : blockwiseDependentIndexPairs) {
        for (const auto &dependentIndexPair : dependentIndexPairs) {
            spdlog::debug("({} {})", dependentIndexPair.first, dependentIndexPair.second);
        }
        spdlog::debug(" ");
    }
    spdlog::debug(" ");


    // Check each environment block
    std::vector<EnvironmentBlock> blocks;
    for (const auto &indexPairsOfBlock : blockwiseDependentIndexPairs) {

        auto block = EnvironmentBlock(indexPairsOfBlock,
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

    // Build Sequence of blocks and test again
    EnvironmentBlockSequence jointBlocks(permutee.molecule_.electrons(), reference.molecule_.electrons());
    for (const auto &block : blocks) {
        auto conservingQ = jointBlocks.addBlock(block, distanceMatrixCovarianceTolerance);

        if (!conservingQ) {
            spdlog::debug("exited early (not inter-block distance preserving)");
            return {earlyExitResult};
        }
    }

    std::vector<BestMatch::DescendingMetricResult> results;

    for (const auto &jointPermutedPermuteeIndices : jointBlocks.jointPermutedPermuteeIndicesCollection_) {
        assert(jointPermutedPermuteeIndices.size() == size_t(N) &&
               "The found index ordering size must match the number of electrons.");

        auto jointReferenceIndices = jointBlocks.jointReferenceIndices_;

        // construct permutations in the kit system
        Eigen::VectorXi finalPermIndicesInKitSystem(N);

        // determine final permutation indices and lowest environmental similarity value for each electron
        double lowestEnvironmentalSimilarityOfParticle = 1.0; // start with maximal value
        for (size_t i = 0; i < N; ++i) {

            auto permIdx = jointPermutedPermuteeIndices[i];
            auto refIdx = jointReferenceIndices[i];

            // create permuation vector
            finalPermIndicesInKitSystem[permIdx] = refIdx;
            spdlog::debug("({} {})", permIdx, refIdx);

            // determine lowest similarity
            auto environmentalSimilarity = environmentalSimilarities(permIdx, refIdx);
            if (environmentalSimilarity < lowestEnvironmentalSimilarityOfParticle)
                lowestEnvironmentalSimilarityOfParticle = environmentalSimilarity;
        }
        spdlog::debug("Permutation (particle-kit system): {}", ToString::vectorXiToString(finalPermIndicesInKitSystem));

        Eigen::PermutationMatrix<Eigen::Dynamic> perm(finalPermIndicesInKitSystem);

        results.emplace_back(BestMatch::DescendingMetricResult(
                {lowestEnvironmentalSimilarityOfParticle, referenceFromKit * perm * permuteeToKit}));
    }
    std::sort(results.begin(), results.end()); // higher metric values come first

    spdlog::set_level(spdlog::level::debug);
    for (auto r : results)
        spdlog::debug("Metric: {}, Permutation: {} (lab system)", r.metric,
                      ToString::vectorXiToString(r.permutation.indices()));
    spdlog::set_level(spdlog::level::info);
    return results;
}


double
BestMatch::SOAPSimilarity::earlyExitMetric(const Eigen::MatrixXd &bestMatchPermutedEnvironmentalSimilarities) {
    //environmentalSimilarities.diagonal().sum() / N; // TODO: average or better smallest component
    return bestMatchPermutedEnvironmentalSimilarities.diagonal().minCoeff();
}

BestMatch::DescendingMetricResult BestMatch::SOAPSimilarity::compare(
        const MolecularSpectrum &permutee,
        const MolecularSpectrum &reference,
        double distanceMatrixCovarianceTolerance, double soapThreshold) {
    return BestMatch::SOAPSimilarity::getBestMatchResults(permutee, reference, distanceMatrixCovarianceTolerance,
                                                          soapThreshold).front();
}

//
//Returns pairs of indices in the kit system for each block of dependent environments from the environment similarity matrix in the kit system
//first index: permutee environment index in the kit system
//second index: reference environment index in the kit system
//
std::vector<std::deque<std::pair<Eigen::Index, Eigen::Index>>>
BestMatch::SOAPSimilarity::getBlockwiseDependentIndexPairs(
        const Eigen::MatrixXd &environmentalSimilarities,
        const Eigen::PermutationMatrix<Eigen::Dynamic> &bestMatch,
        double soapThreshold) {
    std::vector<std::deque<std::pair<Eigen::Index, Eigen::Index>>> listOfDependentIndicesLists;

    auto numberOfEnvironments = environmentalSimilarities.rows();

    auto permutedEnvironments = bestMatch * environmentalSimilarities;// apply environmental best match
    // TODO: giving a reference to the best-match permuted environmental similarity matrix
    //  results in failure of test ListOfDependentIndices4

    std::set<Eigen::Index> indicesToSkip;
    auto epsilon = sqrt(std::numeric_limits<double>::epsilon());

    for (Eigen::Index i = 0; i < numberOfEnvironments; ++i) {

        if (indicesToSkip.find(i) == indicesToSkip.end()) {
            listOfDependentIndicesLists.emplace_back(std::deque({std::pair<Eigen::Index, Eigen::Index>(0, i)}));

            for (Eigen::Index j = i + 1; j < numberOfEnvironments; ++j) {
                if (permutedEnvironments(i, j) >= soapThreshold - epsilon) {
                    indicesToSkip.emplace(j);
                    listOfDependentIndicesLists.back().emplace_back(std::make_pair(0, j));
                }
            }
        }
    }

    // permutee
    Eigen::PermutationMatrix<Eigen::Dynamic> inverseBestMatch = bestMatch.inverse();
    // change indices back to form best match permuted kit system to kit system
    for (auto &list : listOfDependentIndicesLists) {
        auto copy = list;
        for (size_t i = 0; i < list.size(); ++i) {
            list[i].first = inverseBestMatch.indices()[copy[i].second];
        }
    }

    return listOfDependentIndicesLists;
}

/*
 * Calculates the positional distances in the kit system of the electrons specified by the indices
 */
Eigen::MatrixXd
BestMatch::SOAPSimilarity::calculateDistanceCovarianceMatrixOfSelectedIndices(const ElectronsVector &electronsVector,
                                                                              const std::vector<Eigen::Index> &kitSystemIndices) {
    Eigen::MatrixXd positionalDistances(kitSystemIndices.size(), kitSystemIndices.size());

    const auto &positions = electronsVector.positionsVector();

    // the fromKit permutation is needed to find the indices in the Kit system
    auto fromKit = ParticleKit::fromKitPermutation(
            electronsVector).indices(); // TODO refactor into class that stores fromKit and toKit permutations
    assert(fromKit.size() >= kitSystemIndices.size());

    for (Eigen::size_t i = 0; i < kitSystemIndices.size(); ++i)
        for (Eigen::size_t j = 0; j < kitSystemIndices.size(); ++j)
            positionalDistances(i, j) = (positions[fromKit[kitSystemIndices[j]]] -
                                         positions[fromKit[kitSystemIndices[i]]]).norm();

    return positionalDistances;
};
