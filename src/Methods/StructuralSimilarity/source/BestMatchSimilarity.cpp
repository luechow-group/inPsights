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

std::vector<BestMatch::Result> BestMatch::SOAPSimilarity::getBestMatchResults(
        const MolecularSpectrum &permutee,
        const MolecularSpectrum &reference,
        double distanceMatrixCovarianceTolerance, double soapThreshold) {

    assert(ParticleKit::isSubsetQ(permutee.molecule_)
           && "The permutee must be a subset of the particle kit.");
    assert(ParticleKit::isSubsetQ(reference.molecule_)
           && "The reference must be a subset of the particle kit.");

    auto permuteeToKit = ParticleKit::toKitPermutation(permutee.molecule_.electrons());
    auto permuteeFromKit = ParticleKit::toKitPermutation(permutee.molecule_.electrons());
    auto referenceFromKit = ParticleKit::fromKitPermutation(reference.molecule_.electrons());

    // TODO assert that identical number of electrons and same atom geometry? Is this constraint needed? What happens with rows/cols of zero?
    auto N = permutee.molecule_.electrons().numberOfEntities();

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
            ToString::matrixXdToString(environmentalSimilarities,3));

    // find the best match permutation of the environments (of the permutee => rows are permuted) that maximizes the diagonal)
    auto bestMatch = Hungarian<double>::findMatching(environmentalSimilarities, Matchtype::MAX);

    auto bestMatchPermutedEnvironmentalSimilarities = bestMatch*environmentalSimilarities;
    spdlog::debug("Best-Match permutation: {}\n Best-match permuted environmental similarity matrix "
                  "(in the particle kit system, permutee = rows, reference = cols)\n{}",
                  ToString::vectorXiToString(bestMatch.indices()),
            ToString::matrixXdToString(bestMatchPermutedEnvironmentalSimilarities,3));

    // some indices might depend on each other (are equivalent, e.g. two electrons in a nucleus might be swapped)
    // obtain dependent indices in the kit system (blocks of indices might )
    auto blockwiseDependentIndexPairs =
            BestMatch::SOAPSimilarity::getBlockwiseDependentIndexPairs(environmentalSimilarities, bestMatch, soapThreshold);

    std::deque<std::vector<std::deque<std::pair<Eigen::Index,Eigen::Index>>>> intraBlockTestedEnvironmentCombinations;
    for(const auto& blockOfDependentIndexPairs : blockwiseDependentIndexPairs) {
        // vary systematically
        std::vector<std::deque<std::pair<Eigen::Index,Eigen::Index>>> intraBlockTestedEnvironmentCombination;
        BestMatch::SOAPSimilarity::varySimilarEnvironmentsInBlock(permutee.molecule_, reference.molecule_, blockOfDependentIndexPairs, {},
                                                                  intraBlockTestedEnvironmentCombination,
                                                                  distanceMatrixCovarianceTolerance);

        // if one block does not contain a distance preserving permutation, no overall distance preserving permutation can be found
        if(intraBlockTestedEnvironmentCombination.empty()) {
            // no distant preserving permutation could be found
            spdlog::debug("exited early (not intra-block distance preserving");
            auto soapMetric = bestMatchPermutedEnvironmentalSimilarities.diagonal().sum() / N;
            //auto soapMetric = environmentalSimilarities.diagonal().minCoeff() / N; // TODO: average or better smallest component
            return {{soapMetric, bestMatch}};//TODO better return zero or bool?
        } else {
            intraBlockTestedEnvironmentCombinations.emplace_back(intraBlockTestedEnvironmentCombination);
        }
    }

    auto interBlockTestedEnvironmentCombinations = combineBlocks(
            permutee.molecule_, reference.molecule_,
            intraBlockTestedEnvironmentCombinations,
            distanceMatrixCovarianceTolerance);

    if(interBlockTestedEnvironmentCombinations.empty()) {
        // no distant preserving permutation could be found
        auto soapMetric = bestMatchPermutedEnvironmentalSimilarities.diagonal().sum() / N; // TODO WHICH MATRIX SHOULD BE USED
        //auto soapMetric = environmentalSimilarities.diagonal().minCoeff() / N; // TODO: average or better smallest component
        spdlog::debug("exited early (not inter-block distance preserving");
        return {{soapMetric, bestMatch}};//TODO better return zero or bool?
    }

    std::vector<BestMatch::Result> results;

    // Depending on the combination sequence, the metric values can deviate for preserving sequences
    // - thus we have to calculate it for each sequence
    for(const auto& indexPairOfInterBlockCombination : interBlockTestedEnvironmentCombinations) {
        assert(indexPairOfInterBlockCombination.size() == N && "The found index ordering size must match the number of electrons.");

        // construct permutations in the kit system
        Eigen::VectorXi finalPermIndicesInKitSystem(N);

        // determine final permutation
        for (size_t j = 0; j < indexPairOfInterBlockCombination.size(); ++j) {
            spdlog::debug("({} {})", indexPairOfInterBlockCombination[j].first, indexPairOfInterBlockCombination[j].second);
            finalPermIndicesInKitSystem[indexPairOfInterBlockCombination[j].first] = indexPairOfInterBlockCombination[j].second;
        }
        spdlog::debug("Permutation (particle-kit system): {}",ToString::vectorXiToString(finalPermIndicesInKitSystem));

        Eigen::PermutationMatrix<Eigen::Dynamic> perm(finalPermIndicesInKitSystem);

        // find smallest element for metric in the matrix of environmentalSimilarities
        // (NOT in bestMatchPermutedEnvironmentalSimilarities, which only serves the purpose of finding dependent environments)
        assert(!indexPairOfInterBlockCombination.empty());
        double lowestEnvironmentalSimilarityOfParticle = 1.0;

        for(const auto & pair : indexPairOfInterBlockCombination){
            auto environmentalSimilarityOfCurrentParticle = environmentalSimilarities(pair.first, pair.second);
            if(environmentalSimilarityOfCurrentParticle < lowestEnvironmentalSimilarityOfParticle)
                lowestEnvironmentalSimilarityOfParticle = environmentalSimilarityOfCurrentParticle;
        }

        results.emplace_back(BestMatch::Result({lowestEnvironmentalSimilarityOfParticle, referenceFromKit * perm * permuteeToKit}));
    }
    std::sort(results.begin(), results.end());

    for(auto r : results)
        spdlog::debug("Metric: {}, Permutation (lab system): {}", r.metric, ToString::vectorXiToString(r.permutation.indices()));

    return results;
}

BestMatch::Result BestMatch::SOAPSimilarity::compare(
        const MolecularSpectrum &permutee,
        const MolecularSpectrum &reference,
        double distanceMatrixCovarianceTolerance, double soapThreshold) {

    return BestMatch::SOAPSimilarity::getBestMatchResults(permutee, reference, distanceMatrixCovarianceTolerance, soapThreshold).back();
}

/*
 * This method finds combinations of similar environments that are distance preserving across all blocks
 */
std::vector<std::deque<std::pair<Eigen::Index,Eigen::Index>>> BestMatch::SOAPSimilarity::combineBlocks(
        const MolecularGeometry &permutee,
        const MolecularGeometry &reference,
        const std::deque<std::vector<std::deque<std::pair<Eigen::Index,Eigen::Index>>>> &intraBlockDistanceCombinations, /*
        contains blocks of index pairs of dependent similar environments. The first and second index belong to the permutee and reference, respectively.
        For each block, it was checked that the dependent indices are distance preserving under permutation(via the covariance distance matrix difference)*/
        double distanceMatrixCovarianceTolerance){


    assert(!intraBlockDistanceCombinations.empty());

    // This list contains combinations of distance preserving environment combinations of all blocks
    // The first block of intraBlockDistanceCombinations is distance preserving because this was checked initially in varySimilarEnvironments()
    // It grows by adding intraBlockDistanceCombinations that prove to be distance preserving with the existing ones
    std::vector<std::deque<std::pair<Eigen::Index,Eigen::Index>>> interBlockTestedCombinations
    = intraBlockDistanceCombinations.front();

    // all subsequent blocks have to be tested together with the already surviving blocks

    // this loop updates the survivor combination list
    for(auto intraBlockIt = std::next(intraBlockDistanceCombinations.begin());
    intraBlockIt != intraBlockDistanceCombinations.end(); ++intraBlockIt){

        std::vector<std::deque<std::pair<Eigen::Index,Eigen::Index>>> passedInterBlockTestedCombinations;

        // for all combinations of the block that survived till now
        for (auto passingInterBlockTestedCombination = interBlockTestedCombinations.begin();
             passingInterBlockTestedCombination != interBlockTestedCombinations.end(); ++passingInterBlockTestedCombination) {

            // append a ne trial combination from an intra-tested block to a chain of inter-teste combinations
            for (auto &trialIntraBlockTestedCombination : *intraBlockIt) {

                std::deque<Eigen::Index> permuteeIndices, referenceIndices;//(potentialSurvivor.size());
                for(auto i : *passingInterBlockTestedCombination) {
                    permuteeIndices.emplace_back(i.first);
                    referenceIndices.emplace_back(i.second);

                    spdlog::debug("({} {})", i.first, i.second);
                }
                spdlog::debug("_____");
                for(auto i : trialIntraBlockTestedCombination){
                    permuteeIndices.emplace_back(i.first);
                    referenceIndices.emplace_back(i.second);
                    spdlog::debug("({} {})", i.first, i.second);
                }

                // evaluate if the current sequence of combinations is distance preserving
                auto covA = calculateDistanceCovarianceMatrixOfSelectedIndices(permutee.electrons(),
                                                                               permuteeIndices);
                auto covB = calculateDistanceCovarianceMatrixOfSelectedIndices(reference.electrons(),
                                                                               referenceIndices);

                // The maximal distance matrix differences should be smaller than the similarity radius.
                auto covarianceDifference = (covB - covA);
                spdlog::debug("Distance covariance matrix difference:\n{}", ToString::matrixXdToString(covarianceDifference, 3));

                auto conservingQ = covarianceDifference.array().abs().maxCoeff() <= distanceMatrixCovarianceTolerance;
                if (conservingQ) {
                    spdlog::debug(" + distance conserving");
                    auto newlyPassedInterBlockTestedCombination = *passingInterBlockTestedCombination;
                    newlyPassedInterBlockTestedCombination.insert(newlyPassedInterBlockTestedCombination.end(),
                                                                  trialIntraBlockTestedCombination.begin(),
                                                                  trialIntraBlockTestedCombination.end());
                    passedInterBlockTestedCombinations.emplace_back(newlyPassedInterBlockTestedCombination);
                } else {
                    spdlog::debug(" - not distance conserving");
                }
            }
        }
        interBlockTestedCombinations = passedInterBlockTestedCombinations;
    }
    return interBlockTestedCombinations;
}

/*
 * Recursive function!
 * Vary similar environments initially gets a list containing blocks of pairs of dependent indices from the environmental similarity matrix.
 *
 * PURPOSE:
 * Checks the covariance matrix to assert that distances between the selected electrons of a block of similar environments are preserved during their permutation () all
 */
void BestMatch::SOAPSimilarity::varySimilarEnvironmentsInBlock(
        const MolecularGeometry &permutee,
        const MolecularGeometry &reference,
        std::deque<std::pair<Eigen::Index,Eigen::Index>> remainingIndexPairs, // contains a list of dependent indices from a block
        const std::deque<std::pair<Eigen::Index,Eigen::Index>>& survivingIndexPairs, // is empty in the beginning
        std::vector<std::deque<std::pair<Eigen::Index,Eigen::Index>>> &distancePreservingEnvironmentCombinations,
        double distanceMatrixCovarianceTolerance) {

    if(remainingIndexPairs.empty()) {
        distancePreservingEnvironmentCombinations.emplace_back(survivingIndexPairs);
        return;
    }

    // iterate over blocks of dependent index pairs
    for (auto it = remainingIndexPairs.begin(); it != remainingIndexPairs.end(); ++it){
        auto potentialSurvivingIndexPairs = survivingIndexPairs; // this gets longer during recursion or stops

        potentialSurvivingIndexPairs.emplace_back(*it);

        // collect all indices that are dependent in the block for both, the permutee and the reference
        std::deque<Eigen::Index> potentiallySurvingPermuteeIndices, potentiallySurvivingReferenceIndices;
        for(auto i : potentialSurvivingIndexPairs) {
            potentiallySurvingPermuteeIndices.emplace_back(i.first);
            potentiallySurvivingReferenceIndices.emplace_back(i.second);
        }
        //TODO is sorting necessary? //TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
        // TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
        // CHECK FOR OTHER TESTS
        //std::sort(potentiallySurvivingReferenceIndices.begin(), potentiallySurvivingReferenceIndices.end());

        // DO WE LOOSE SOME COMBINATIONS HERE BY JUST TRYING ONE COMBINATION OF THE BLOCK?
        // WHERE IS THE VARIATION

        auto covA = calculateDistanceCovarianceMatrixOfSelectedIndices(permutee.electrons(), potentiallySurvingPermuteeIndices);
        auto covB = calculateDistanceCovarianceMatrixOfSelectedIndices(reference.electrons(), potentiallySurvivingReferenceIndices);

        // check if distance covariance is preserved for all combinations of the dependent indices //TODO do we potentially loose conbinations here?
        auto conservingQ = (covB-covA).array().abs().maxCoeff() <= distanceMatrixCovarianceTolerance;

        if(conservingQ) { // go deeper
            auto updatedRemainingIndexPairsOfAllBlocks = remainingIndexPairs;

            updatedRemainingIndexPairsOfAllBlocks.erase(updatedRemainingIndexPairsOfAllBlocks.begin() + std::distance(remainingIndexPairs.begin(), it));

            varySimilarEnvironmentsInBlock(permutee, reference, updatedRemainingIndexPairsOfAllBlocks,
                                           potentialSurvivingIndexPairs,
                                           distancePreservingEnvironmentCombinations, distanceMatrixCovarianceTolerance);
        }
    }
};

/*
 * Returns pairs of indices in the kit system for each block of dependent environments from the environment similarity matrix in the kit system
 * first index: permutee environment index in the kit system
 * second index: reference environment index in the kit system
 */
std::vector<std::deque<std::pair<Eigen::Index,Eigen::Index>>>
BestMatch::SOAPSimilarity::getBlockwiseDependentIndexPairs(
        const Eigen::MatrixXd &environmentalSimilarities,
        const Eigen::PermutationMatrix<Eigen::Dynamic> &bestMatch,
        double soapThreshold) {
    std::vector<std::deque<std::pair<Eigen::Index,Eigen::Index>>> listOfDependentIndicesLists;

    auto numberOfEnvironments = environmentalSimilarities.rows();

    auto permutedEnvironments = bestMatch*environmentalSimilarities;// apply environmental best match

    std::set<Eigen::Index> indicesToSkip;
    auto epsilon = sqrt(std::numeric_limits<double>::epsilon());

    for (Eigen::Index i = 0; i < numberOfEnvironments; ++i) {

        if (indicesToSkip.find(i) == indicesToSkip.end()) {
            std::deque temp = {std::pair<Eigen::Index, Eigen::Index>(0,i)};
            listOfDependentIndicesLists.emplace_back(temp); // TODO refactor: remove temp

            for (Eigen::Index j = i + 1; j < numberOfEnvironments; ++j) {
                if (permutedEnvironments(i, j) >= soapThreshold - epsilon) {
                    indicesToSkip.emplace(j);
                    listOfDependentIndicesLists.back().emplace_back(std::make_pair(0,j));
                }
            }
        }
    }

    // permutee
    Eigen::PermutationMatrix<Eigen::Dynamic> inverseBestMatch = bestMatch.inverse();
    // change indices back to form best match permuted kit system to kit system
    for (auto & list : listOfDependentIndicesLists) {
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
Eigen::MatrixXd BestMatch::SOAPSimilarity::calculateDistanceCovarianceMatrixOfSelectedIndices(const ElectronsVector &electronsVector,
                                                                                              std::deque<Eigen::Index> kitSystemIndices){
    Eigen::MatrixXd positionalDistances(kitSystemIndices.size(), kitSystemIndices.size());

    const auto& positions = electronsVector.positionsVector();

    // the fromKit permutation is needed to find the indices in the Kit system
    auto fromKit = ParticleKit::fromKitPermutation(electronsVector).indices(); // TODO refactor into class that stores fromKit and toKit permutations

    assert(fromKit.size() >= kitSystemIndices.size());

    for (Eigen::size_t i = 0; i < kitSystemIndices.size(); ++i)
        for (Eigen::size_t j = 0; j < kitSystemIndices.size(); ++j)
            positionalDistances(i, j) = (positions[fromKit[kitSystemIndices[j]]] - positions[fromKit[kitSystemIndices[i]]]).norm();

    return positionalDistances;
};
