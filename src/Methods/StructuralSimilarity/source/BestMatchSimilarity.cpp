//
// Created by heuer on 03.04.19.
//

#include "BestMatchSimilarity.h"
#include <LocalSimilarity.h>
#include <Hungarian.h>
#include <BestMatchSimilarity.h>
#include <Metrics.h>
#include <Combinatorics.h>
#include <limits>
#include <deque>
#include <set>
#include <vector>
#include <numeric> //std::iota

#include <iomanip>

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

std::vector<BestMatch::Result> BestMatch::SOAPSimilarity::getAllBestMatchResults(
        const MolecularSpectrum &permutee,
        const MolecularSpectrum &reference,
        double similarityRadius, double soapThreshold) {

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

    Eigen::MatrixXd environmentalSimilarities = calculateEnvironmentalSimilarityMatrix(permutee, reference);

    Eigen::PermutationMatrix<Eigen::Dynamic> bestMatch
            = Hungarian<double>::findMatching(environmentalSimilarities, Matchtype::MAX);

    auto dependentIndicesLists =
            BestMatch::SOAPSimilarity::getListOfDependentIndicesLists(environmentalSimilarities, soapThreshold);

    std::deque<std::vector<std::deque<Eigen::Index>>> distancePreservingEnvironmentCombinationsOfAllBlocks;

    for(const auto& list : dependentIndicesLists) {
        std::vector<std::deque<Eigen::Index>> distancePreservingEnvironmentCombinations;
        BestMatch::SOAPSimilarity::varySimilarEnvironments(permutee.molecule_, reference.molecule_,list, {},
                                                           distancePreservingEnvironmentCombinations, similarityRadius);

        if(distancePreservingEnvironmentCombinations.empty()) {
            // no distant preserving permutation could be found f
            auto soapMetric = environmentalSimilarities.diagonal().sum() / N;
            return {{soapMetric, bestMatch}};
        }
        distancePreservingEnvironmentCombinationsOfAllBlocks.emplace_back(distancePreservingEnvironmentCombinations);
    }

    auto survivingDistancePreservingEnvironmentCombinations = combineBlocks(
            permutee.molecule_, reference.molecule_,
            distancePreservingEnvironmentCombinationsOfAllBlocks,
            similarityRadius);

    if(survivingDistancePreservingEnvironmentCombinations.empty()) {
        // no distant preserving permutation could be found
        auto soapMetric = environmentalSimilarities.diagonal().sum() / N;
        return {{soapMetric, bestMatch}};
    }

    auto indexReorderingPermutation
            = obtainIndexReorderingPermutationOverAllBlocks(distancePreservingEnvironmentCombinationsOfAllBlocks);

    // bring the surviving distance preserving environment combinations back to the original order
    for(auto& indices : survivingDistancePreservingEnvironmentCombinations) {
        auto copy = indices;

        for (std::size_t i = 0; i < indices.size(); ++i) {
            indices[i] = copy[indexReorderingPermutation[i]];
        }
    }

    std::vector<BestMatch::Result> results;

    for(auto& foundIndicesOrder : survivingDistancePreservingEnvironmentCombinations) {
        assert(foundIndicesOrder.size() == N && "The found index ordering size must match the number of electrons.");

        Eigen::VectorXi permIndices(N);
        for (std::size_t i = 0; i < foundIndicesOrder.size(); ++i) {
            permIndices[i] = foundIndicesOrder[i];
        }
        Eigen::PermutationMatrix<Eigen::Dynamic> perm(permIndices);

        // find smallest element
        assert(foundIndicesOrder.size() > 0);
        double smallestEnvironmentalBestMatchSimilarity = environmentalSimilarities(foundIndicesOrder[0],0);
        for (std::size_t i = 1; i < foundIndicesOrder.size(); ++i) {
            auto ithEnvironmentalBestMatchSimilarity = environmentalSimilarities(foundIndicesOrder[i],i);
            if(ithEnvironmentalBestMatchSimilarity < smallestEnvironmentalBestMatchSimilarity)
                smallestEnvironmentalBestMatchSimilarity = ithEnvironmentalBestMatchSimilarity;
        }

        results.emplace_back(BestMatch::Result({smallestEnvironmentalBestMatchSimilarity, referenceFromKit * perm * permuteeToKit}));
    }
    std::sort(results.begin(), results.end());

    return results;
}

BestMatch::Result BestMatch::SOAPSimilarity::compare(
        const MolecularSpectrum &permutee,
        const MolecularSpectrum &reference,
        double similarityRadius, double soapThreshold) {

    return BestMatch::SOAPSimilarity::getAllBestMatchResults(permutee,reference,similarityRadius, soapThreshold).back();
}

std::vector<std::deque<Eigen::Index>> BestMatch::SOAPSimilarity::combineBlocks(
        const MolecularGeometry &permutee,
        const MolecularGeometry &reference,
        const std::deque<std::vector<std::deque<Eigen::Index>>> &distancePreservingEnvironmentCombinationsOfRemainingBlocks,
        double similarityRadius){

    assert(!distancePreservingEnvironmentCombinationsOfRemainingBlocks.empty());

    std::vector<std::deque<Eigen::Index>> survivors = distancePreservingEnvironmentCombinationsOfRemainingBlocks.front();

    for(auto blockIt = std::next(distancePreservingEnvironmentCombinationsOfRemainingBlocks.begin()); blockIt != distancePreservingEnvironmentCombinationsOfRemainingBlocks.end(); ++blockIt){

        std::vector<std::deque<Eigen::Index>> newSurvivors;

        // for every remaining distance preserving environment combinations
        for (auto potentialSurvivorIt = survivors.begin(); potentialSurvivorIt != survivors.end(); ++potentialSurvivorIt) {

            // check all combinations in the current block
            for (auto &distancePreservingEnvironmentCombination : *blockIt) {

                // TRICK
                // the potential survivor is already distance preserving, now add the new candidates original index order
                // appending the new indices is fine since they are independent from the ones before


                auto newPotentialSurvivor = *potentialSurvivorIt;
                auto originalIndexOrder = *potentialSurvivorIt;

                // append the new distance preserving combination of the current block
                newPotentialSurvivor.insert(newPotentialSurvivor.end(),
                                              distancePreservingEnvironmentCombination.begin(),
                                              distancePreservingEnvironmentCombination.end());
                //hacky: this uses the fact, that the original index order is ascending


                // compare the potential survivor with the unpermuted reference.
                // The distance matrix differences should be zero.
                auto newCandidateOriginalOrder = distancePreservingEnvironmentCombination;

                originalIndexOrder.insert(originalIndexOrder.end(),
                                          newCandidateOriginalOrder.begin(),
                                          newCandidateOriginalOrder.end());
                std::sort(originalIndexOrder.begin(),
                          originalIndexOrder.end());

                auto covA = indicesBlockCovariance(permutee.electrons(), newPotentialSurvivor);
                auto covB = indicesBlockCovariance(reference.electrons(), originalIndexOrder);

                auto conservingQ = (covB - covA).array().abs().maxCoeff() <= similarityRadius;

                if (conservingQ)
                    newSurvivors.emplace_back(newPotentialSurvivor);
            }
        }
        survivors = newSurvivors;
    }
    return survivors;
}

std::vector<Eigen::Index> BestMatch::SOAPSimilarity::obtainIndexReorderingPermutationOverAllBlocks(
        const std::deque<std::vector<std::deque<Eigen::Index>>> &distancePreservingEnvironmentCombinationsOfRemainingBlocks) {

    // extract blockwise ordered indices
    std::vector<Eigen::Index> blockWiseReorderedIndices;
    for (std::size_t i = 0; i < distancePreservingEnvironmentCombinationsOfRemainingBlocks.size(); ++i) {

        auto examplaricOrder = distancePreservingEnvironmentCombinationsOfRemainingBlocks[i].front();
        std::sort(examplaricOrder.begin(), examplaricOrder.end());

        for (auto index : examplaricOrder) {
            blockWiseReorderedIndices.emplace_back(index);
        }
    }

    // initialize identity permutation
    std::vector<Eigen::Index> overallReorderedIndices(blockWiseReorderedIndices.size());
    std::iota(overallReorderedIndices.begin(), overallReorderedIndices.end(), 0);

    // sort and store index permutation
    sort(overallReorderedIndices.begin(), overallReorderedIndices.end(),
         [&](const int &a, const int &b) {
             return (blockWiseReorderedIndices[a] < blockWiseReorderedIndices[b]);
         }
    );

    return overallReorderedIndices;
};


void BestMatch::SOAPSimilarity::varySimilarEnvironments(
        const MolecularGeometry &permutee,
        const MolecularGeometry &reference,
        std::deque<Eigen::Index> dependentIndices, // contains list of dependent indices
        std::deque<Eigen::Index> surviving, // is empty in the beginning
        std::vector<std::deque<Eigen::Index>> &allPerms,
        double similarityRadius) {


    if(dependentIndices.size() == 0) {
        //surviving.emplace_back(dependentIndices.front());
        allPerms.emplace_back(surviving);
        return;
    }

    for (auto it = dependentIndices.begin(); it != dependentIndices.end(); ++it){
        auto potentialSurvivor = surviving;

        potentialSurvivor.emplace_back(*it);

        auto covA = indicesBlockCovariance(permutee.electrons(), potentialSurvivor);

        auto originalIndexOrder = potentialSurvivor;
        std::sort(originalIndexOrder.begin(),originalIndexOrder.end()); //TODO hacky: this uses the fact, that the original index order is ascending

        auto covB = indicesBlockCovariance(reference.electrons(), originalIndexOrder);

        auto conservingQ = (covB-covA).array().abs().maxCoeff() <= similarityRadius;

        if(conservingQ) { // go deeper
            auto remainingCopy = dependentIndices;

            remainingCopy.erase(remainingCopy.begin() + std::distance(dependentIndices.begin(), it));

            varySimilarEnvironments(permutee, reference, remainingCopy, potentialSurvivor, allPerms, similarityRadius);
        }
    }
};

std::vector<std::deque<Eigen::Index>>
BestMatch::SOAPSimilarity::getListOfDependentIndicesLists(
        const Eigen::MatrixXd &environmentalSimilarities, double soapThreshold) {
    std::vector<std::deque<Eigen::Index>> listOfDependentIndicesLists;

    auto numberOfEnvironments = environmentalSimilarities.rows();

    std::set<Eigen::Index> indicesToSkip;
    auto epsilon = sqrt(std::numeric_limits<double>::epsilon());

    for (Eigen::Index i = 0; i < numberOfEnvironments; ++i) {

        if (indicesToSkip.find(i) == indicesToSkip.end()) {
            listOfDependentIndicesLists.emplace_back(std::deque({i}));

            for (Eigen::Index j = i + 1; j < numberOfEnvironments; ++j) {
                if (environmentalSimilarities(i, j) >= soapThreshold - epsilon) {
                    indicesToSkip.emplace(j);
                    listOfDependentIndicesLists.back().emplace_back(j);
                }
            }
        }
    }
    return listOfDependentIndicesLists;
}

Eigen::MatrixXd BestMatch::SOAPSimilarity::indicesBlockCovariance(const ElectronsVector &electronsVector,
                                                                  std::deque<Eigen::Index> indices){
    Eigen::MatrixXd distDiffBlock(indices.size(), indices.size());

    const auto& positions = electronsVector.positionsVector();

    // the toKit permutation is needed to find the indices in the Kit system
    auto permIndices = ParticleKit::toKitPermutation(electronsVector).indices();

    for (Eigen::size_t i = 0; i < indices.size(); ++i)
        for (Eigen::size_t j = 0; j < indices.size(); ++j)
            distDiffBlock(i,j) = (positions[permIndices[indices[j]]] - positions[permIndices[indices[i]]]).norm();

    return distDiffBlock;
};
