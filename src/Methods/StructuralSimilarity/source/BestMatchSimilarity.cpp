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
#include <Eigen/Core>
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

std::vector<BestMatch::Result> BestMatch::SOAPSimilarity::getBestMatchResults(
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
    std::cout<< std::endl << std::setprecision(5) << environmentalSimilarities << std::endl;

    auto bestMatch
    = Hungarian<double>::findMatching(environmentalSimilarities, Matchtype::MAX);

    std::cout<< std::endl << std::setprecision(5) << bestMatch.indices().transpose() << std::endl;

    std::cout<< std::endl << std::setprecision(5) <<ParticleKit::toKitPermutation(permutee.molecule_.electrons()).indices().transpose() << std::endl;
    std::cout<< std::endl << std::setprecision(5) <<ParticleKit::toKitPermutation(reference.molecule_.electrons()).indices().transpose() << std::endl;

    std::cout<< std::endl << std::setprecision(5) << bestMatch*environmentalSimilarities << std::endl;

    // obtain indices in the kit system
    auto dependentIndicesLists =
            BestMatch::SOAPSimilarity::getListOfDependentIndicesLists(environmentalSimilarities, bestMatch, soapThreshold);

    for(auto i : dependentIndicesLists) {
        for (auto j : i) {
            std::cout << "("<< j.first << "," << j.second << ") ";
        }
        std::cout << std::endl;
    }

    std::deque<std::vector<std::deque<std::pair<Eigen::Index,Eigen::Index>>>> distancePreservingEnvironmentCombinationsOfAllBlocks;

    for(const auto& list : dependentIndicesLists) {
        // vary systematically
        std::vector<std::deque<std::pair<Eigen::Index,Eigen::Index>>> distancePreservingEnvironmentCombinations;
        BestMatch::SOAPSimilarity::varySimilarEnvironments(permutee.molecule_, reference.molecule_,list, {},
                                                           distancePreservingEnvironmentCombinations, similarityRadius);

        if(distancePreservingEnvironmentCombinations.empty()) {
            // no distant preserving permutation could be found
            std::wcerr<< "exited early 1" << std::endl;
            auto soapMetric = environmentalSimilarities.diagonal().sum() / N;
            return {{soapMetric, bestMatch}};//TODO better return zero or bool?
        }
        distancePreservingEnvironmentCombinationsOfAllBlocks.emplace_back(distancePreservingEnvironmentCombinations);
    }


    for(auto i : distancePreservingEnvironmentCombinationsOfAllBlocks) {
        for (auto j : i) {
            for (auto k : j) {
                std::cout << "(" << k.first << "," << k.second << ") ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

    auto survivingDistancePreservingEnvironmentCombinations = combineBlocks(
            permutee.molecule_, reference.molecule_,
            distancePreservingEnvironmentCombinationsOfAllBlocks,
            similarityRadius);

    if(survivingDistancePreservingEnvironmentCombinations.empty()) {
        // no distant preserving permutation could be found
        auto soapMetric = environmentalSimilarities.diagonal().sum() / N;
        std::wcerr<< "exited early 2" << std::endl;
        return {{soapMetric, bestMatch}};//TODO better return zero or bool?
    }


    // Reverse the blockwise ordering of permutations
    auto indexReorderingPermutation
            = obtainIndexReorderingPermutationOverAllBlocks(distancePreservingEnvironmentCombinationsOfAllBlocks);

    // bring the surviving distance preserving environment combinations back to the original order
    for(auto& indices : survivingDistancePreservingEnvironmentCombinations) {
        auto copy = indices;
        for (std::size_t i = 0; i < indices.size(); ++i) {
            indices[i] = copy[indexReorderingPermutation[i]];
        }
    }
    std::wcerr<< "index reordering :";
    for (std::size_t i = 0; i < indexReorderingPermutation.size(); ++i) {
        std::cout << indexReorderingPermutation[i] << " ";
    }
    std::cout << std::endl;


    std::vector<BestMatch::Result> results;

    for(auto& foundIndicesOrder : survivingDistancePreservingEnvironmentCombinations) {
        assert(foundIndicesOrder.size() == N && "The found index ordering size must match the number of electrons.");

        Eigen::VectorXi permIndices(N);
        for (std::size_t i = 0; i < foundIndicesOrder.size(); ++i) {
            std::cout << foundIndicesOrder[i].first << " ";
            permIndices[i] = foundIndicesOrder[i].first;
        }
        std::cout << std::endl;
        Eigen::PermutationMatrix<Eigen::Dynamic> perm(permIndices);

        // find smallest element
        assert(foundIndicesOrder.size() > 0);
        double smallestEnvironmentalBestMatchSimilarity = environmentalSimilarities(foundIndicesOrder[0].first,0);
        for (std::size_t i = 1; i < foundIndicesOrder.size(); ++i) {
            auto ithEnvironmentalBestMatchSimilarity = environmentalSimilarities(foundIndicesOrder[i].first,i);
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

    return BestMatch::SOAPSimilarity::getBestMatchResults(permutee, reference, similarityRadius, soapThreshold).back();
}

std::vector<std::deque<std::pair<Eigen::Index,Eigen::Index>>> BestMatch::SOAPSimilarity::combineBlocks(
        const MolecularGeometry &permutee,
        const MolecularGeometry &reference,
        const std::deque<std::vector<std::deque<std::pair<Eigen::Index,Eigen::Index>>>> &distancePreservingEnvironmentCombinationsOfRemainingBlocks,
        double similarityRadius){

    assert(!distancePreservingEnvironmentCombinationsOfRemainingBlocks.empty());

    std::vector<std::deque<std::pair<Eigen::Index,Eigen::Index>>> survivors = distancePreservingEnvironmentCombinationsOfRemainingBlocks.front();

    for(auto blockIt = std::next(distancePreservingEnvironmentCombinationsOfRemainingBlocks.begin()); blockIt != distancePreservingEnvironmentCombinationsOfRemainingBlocks.end(); ++blockIt){

        std::vector<std::deque<std::pair<Eigen::Index,Eigen::Index>>> newSurvivors;

        // for every remaining distance preserving environment combinations
        for (auto potentialSurvivorIt = survivors.begin(); potentialSurvivorIt != survivors.end(); ++potentialSurvivorIt) {

            // check all combinations in the current block
            for (auto &distancePreservingEnvironmentCombination : *blockIt) {

                // TRICK
                // the potential survivor is already distance preserving, now add the new candidates original index order
                // appending the new indices is fine since they are independent from the ones before


                /*auto newPotentialSurvivor = *potentialSurvivorIt;
                auto originalIndexOrder = *potentialSurvivorIt;

                // append the new distance preserving combination of the current block
                newPotentialSurvivor.insert(newPotentialSurvivor.end(),
                                              distancePreservingEnvironmentCombination.begin(),
                                              distancePreservingEnvironmentCombination.end());
                //hacky: this uses the fact, that the original index order is ascending


                // compare the potential survivor with the unpermuted reference.
                originalIndexOrder.insert(originalIndexOrder.end(),
                                          distancePreservingEnvironmentCombination.begin(),
                                          distancePreservingEnvironmentCombination.end());
                std::sort(originalIndexOrder.begin(),
                          originalIndexOrder.end());*/


                std::deque<Eigen::Index> permuteeIndices, referenceIndices;//(potentialSurvivor.size());
                for(auto i : *potentialSurvivorIt) {
                    permuteeIndices.emplace_back(i.first);
                    referenceIndices.emplace_back(i.second);
                }
                for(auto i : distancePreservingEnvironmentCombination){
                    permuteeIndices.emplace_back(i.first);
                    referenceIndices.emplace_back(i.second);
                }
                std::sort(referenceIndices.begin(),referenceIndices.end()); // necessary?
                std::cout << " permutee: ";
                for(auto i : permuteeIndices)
                    std::cout << i << " ";
                std::cout << std::endl << "reference: ";
                for(auto i : referenceIndices)
                    std::cout << i << " ";
                std::cout << std::endl;

                auto covA = indicesBlockCovariance(permutee.electrons(), permuteeIndices);
                auto covB = indicesBlockCovariance(reference.electrons(), referenceIndices);

                // The maximal distance matrix differences should be smaller than the similarity radius.
                auto conservingQ = (covB - covA).array().abs().maxCoeff() <= similarityRadius;
                std::cout << (covB - covA).array().abs() << std::endl;
                if (conservingQ) {
                    auto newSurvivor = *potentialSurvivorIt;
                    newSurvivor.insert(newSurvivor.end(),
                                                distancePreservingEnvironmentCombination.begin(),
                                                distancePreservingEnvironmentCombination.end());
                    newSurvivors.emplace_back(newSurvivor);
                }
            }
        }
        survivors = newSurvivors;
    }
    return survivors;
}

std::vector<Eigen::Index> BestMatch::SOAPSimilarity::obtainIndexReorderingPermutationOverAllBlocks(
        const std::deque<std::vector<std::deque<std::pair<Eigen::Index,Eigen::Index>>>> &distancePreservingEnvironmentCombinationsOfRemainingBlocks) {

    // extract blockwise ordered indices
    std::vector<std::pair<Eigen::Index,Eigen::Index>> blockWiseReorderedIndices;
    for (std::size_t i = 0; i < distancePreservingEnvironmentCombinationsOfRemainingBlocks.size(); ++i) {

        // pick any distance preserving permutation
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
         [&](const Eigen::Index &a, const Eigen::Index &b) {
             return (blockWiseReorderedIndices[a].first < blockWiseReorderedIndices[b].second);
         }
    );

    return overallReorderedIndices;
};


void BestMatch::SOAPSimilarity::varySimilarEnvironments(
        const MolecularGeometry &permutee,
        const MolecularGeometry &reference,
        std::deque<std::pair<Eigen::Index,Eigen::Index>> remaining, // contains list of dependent indices
        std::deque<std::pair<Eigen::Index,Eigen::Index>> surviving, // is empty in the beginning
        std::vector<std::deque<std::pair<Eigen::Index,Eigen::Index>>> &allPerms,
        double similarityRadius) {

    if(remaining.empty()) {
        allPerms.emplace_back(surviving);
        std::cout << "done" << std::endl;
        return;
    }

    for (auto it = remaining.begin(); it != remaining.end(); ++it){
        auto potentialSurvivor = surviving;

        std::cout <<std::endl << "surviving ";
        for(auto i : surviving)
            std::cout << "(" << i.first << "," << i.second << ") ";
        std::cout <<std::endl << "remaining ";

        for(auto i : remaining)
            std::cout << "(" << i.first << "," << i.second << ") ";
        std::cout <<std::endl;

        potentialSurvivor.emplace_back(*it);
        std::cout << "potential survivor ";
        for(auto i : potentialSurvivor)
            std::cout << "(" << i.first << "," << i.second << ") ";
        std::cout <<std::endl;

        //TODO create indices list
        // permutee

        std::deque<Eigen::Index> permuteeIndices, referenceIndices;
        for(auto i : potentialSurvivor) {
            permuteeIndices.emplace_back(i.first);
            referenceIndices.emplace_back(i.second);
        }
        //TODO is sorting necessary?
        std::sort(referenceIndices.begin(),referenceIndices.end());

        std::cout << "permutee " << std::endl;
        auto covA = indicesBlockCovariance(permutee.electrons(), permuteeIndices);
        std::cout << "reference" << std::endl;
        auto covB = indicesBlockCovariance(reference.electrons(), referenceIndices);

        std::cout << (covB-covA).array().abs()<<std::endl;

        auto conservingQ = (covB-covA).array().abs().maxCoeff() <= similarityRadius;

        if(conservingQ) { // go deeper
            auto remainingCopy = remaining;

            remainingCopy.erase(remainingCopy.begin() + std::distance(remaining.begin(), it));

            varySimilarEnvironments(permutee, reference, remainingCopy, potentialSurvivor, allPerms, similarityRadius);
        }
    }
};

std::vector<std::deque<std::pair<Eigen::Index,Eigen::Index>>>
BestMatch::SOAPSimilarity::getListOfDependentIndicesLists(
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
            listOfDependentIndicesLists.emplace_back(temp);


            for (Eigen::Index j = i + 1; j < numberOfEnvironments; ++j) {
                if (permutedEnvironments(i, j) >= soapThreshold - epsilon) {
                    indicesToSkip.emplace(j);
                    listOfDependentIndicesLists.back().emplace_back(std::make_pair(0,j));
                }
            }
        }
    }
    
    // Achtung alles ist im kit system
    
    // Beim variieren der ähnlichen umgebungen des permutee werden die zugehörigen umgebungen der reference benötigt.
    /*
     * Bsp: (0,2) im permutee muss nicht der umgebung (0,2) in der Reference entsprechen -> in diesem fall würden keine matches gefunden
     * Lsg: zusammengehörige umgebungen speichern. Was passiert bei mehreren möglichkeiten?
     * die ref. soll nicht verändert werden - beim blockweisen vergleich soll aber k
     */

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

Eigen::MatrixXd BestMatch::SOAPSimilarity::indicesBlockCovariance(const ElectronsVector &electronsVector,
                                                                  std::deque<Eigen::Index> indices){
    Eigen::MatrixXd positionalDistances(indices.size(), indices.size());

    const auto& positions = electronsVector.positionsVector();

    // the fromKit permutation is needed to find the indices in the Kit system
    auto fromKit = ParticleKit::fromKitPermutation(electronsVector).indices();

    assert(fromKit.size() >= indices.size());

    std::cout << "  kit indices: ";
    for (Eigen::size_t i = 0; i < indices.size(); ++i) {
        std::cout << indices[i] << " ";
    }
    std::cout << std::endl;
    std::cout << "  lab indices: ";
    for (Eigen::size_t i = 0; i < indices.size(); ++i) {
        std::cout << fromKit[indices[i]]<< " ";
    }
    std::cout << std::endl;

    for (Eigen::size_t i = 0; i < indices.size(); ++i)
        for (Eigen::size_t j = 0; j < indices.size(); ++j)
            positionalDistances(i, j) = (positions[fromKit[indices[j]]] - positions[fromKit[indices[i]]]).norm();

    return positionalDistances;
};
