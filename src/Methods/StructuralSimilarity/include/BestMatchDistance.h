//
// Created by heuer on 03.04.19.
//

#ifndef INPSIGHTS_BESTMATCHDISTANCE_H
#define INPSIGHTS_BESTMATCHDISTANCE_H

#include "Hungarian.h"
#include <BestMatch.h>
#include "BestMatchDistance.h"
#include <Metrics.h>
#include <ParticlesVector.h>
#include <Interval.h>

namespace BestMatch {
    namespace Distance {
        template<int positionalNorm = 2>
        Eigen::PermutationMatrix<Eigen::Dynamic> findPermutation(
                const PositionsVector &permutee,
                const PositionsVector &reference) {
            assert(permutee.numberOfEntities() == reference.numberOfEntities());

            auto costMatrix = Metrics::positionalDistances<positionalNorm>(permutee, reference);

            return Hungarian<double>::findMatching(costMatrix);
        }

        template <typename Type>
        Eigen::PermutationMatrix<Eigen::Dynamic> findTypeSeparatingPermutation(
                const ParticlesVector<Type>& particlesVector){

            Eigen::VectorXi typeSerparatingPermutationIndices(particlesVector.numberOfEntities());

            Eigen::Index i = 0;

            for(const auto& [type, count] : particlesVector.typesVector().countTypes()) {

                for (std::size_t j = 0; j < count; ++j) {
                    auto [foundQ, index] = particlesVector.typesVector().findIndexOfEnumeratedType(
                            EnumeratedType<Type>(type,j));
                    assert(foundQ);
                    typeSerparatingPermutationIndices[i] = index;
                    ++i;
                }
            }
            return Eigen::PermutationMatrix<Eigen::Dynamic>(typeSerparatingPermutationIndices);
        }

        template<typename Type, int positionalNorm = 2>
        Eigen::PermutationMatrix<Eigen::Dynamic>
        findTypeSpecificPermutation(ParticlesVector<Type> permutee, ParticlesVector<Type> reference) {
            //TODO//assert(permutee.typesVector() == reference.typesVector()
            //       && "The type vectors must be identical as we assume ordered alpha and beta electrons.");

            assert(permutee.numberOfEntities() == reference.numberOfEntities()
                   && "The number of particles must be identical.");

            auto countedTypes = permutee.typesVector().countTypes();

            auto permuteeSeparatingPermutation = findTypeSeparatingPermutation<Type>(permutee);
            auto referenceSeparatingPermutation = findTypeSeparatingPermutation<Type>(reference);


            Eigen::PermutationMatrix<Eigen::Dynamic> combinedPerm(0);
            Eigen::Index accumulatedCounts = 0;

            for(const auto &typeCount : countedTypes) {

                // Construct cost matrices
                std::list<long> selectedPermuteeIndices,selectedReferenceIndices;

                for (Eigen::Index i = 0; i < typeCount.number_; ++i) {
                    selectedPermuteeIndices.emplace_back(permuteeSeparatingPermutation.indices()[accumulatedCounts+i]);
                    selectedReferenceIndices.emplace_back(referenceSeparatingPermutation.indices()[accumulatedCounts+i]);
                }

                auto costMatrix = Metrics::positionalDistances<positionalNorm>(
                        permutee[selectedPermuteeIndices].positionsVector(),
                        reference[selectedReferenceIndices].positionsVector());

                // find best match and combine it
                combinedPerm = combinePermutations(combinedPerm, Hungarian<double>::findMatching(costMatrix));

                accumulatedCounts += typeCount.number_;
            }

            return referenceSeparatingPermutation * combinedPerm * permuteeSeparatingPermutation.inverse();
        };

        template<int positionalNorm = 2>
        Eigen::PermutationMatrix<Eigen::Dynamic>
        findSpinSpecificPermutation(const ElectronsVector &permutee, const ElectronsVector &reference, bool flipSpinsQ = false) {
                assert(permutee.typesVector() == reference.typesVector()
                       && "The typesvectors must be identical as we assume ordered alpha and beta electrons.");
                assert(permutee.positionsVector().numberOfEntities() == reference.positionsVector().numberOfEntities()
                       && "The number of positions must be identical.");

                auto nAlpha = permutee.typesVector().countOccurence(Spin::alpha);
                auto nBeta = permutee.typesVector().countOccurence(Spin::beta);

                Interval permuteeAlpha{0, nAlpha}, referenceAlpha{};
                Interval permuteeBeta{nAlpha, nBeta}, referenceBeta{};

                if (!flipSpinsQ) {
                    referenceAlpha = permuteeAlpha;
                    referenceBeta = permuteeBeta;
                } else {
                    assert(nAlpha == nBeta
                           && "The number of alpha and beta electrons must match to allow for spin flips.");

                    // flip slice intervals
                    referenceAlpha = permuteeBeta;
                    referenceBeta = permuteeAlpha;
                };

                auto costMatrixAlpha = Metrics::positionalDistances<positionalNorm>(
                        PositionsVector(permutee.positionsVector().asEigenVector().segment(
                                        permuteeAlpha.start() * 3,permuteeAlpha.numberOfEntities() * 3)),
                        PositionsVector(reference.positionsVector().asEigenVector().segment(
                                        referenceAlpha.start() * 3, referenceAlpha.numberOfEntities() * 3)));

                auto bestMatchAlpha = Hungarian<double>::findMatching(costMatrixAlpha);

                auto costMatrixBeta = Metrics::positionalDistances<positionalNorm>(
                        PositionsVector(permutee.positionsVector().asEigenVector().segment(
                                        permuteeBeta.start() * 3, permuteeBeta.numberOfEntities() * 3)),
                        PositionsVector(reference.positionsVector().asEigenVector().segment(
                                        referenceBeta.start() * 3, referenceBeta.numberOfEntities() * 3)));

                auto bestMatchBeta = Hungarian<double>::findMatching(costMatrixBeta);

                if (!flipSpinsQ)
                    return combinePermutations(bestMatchAlpha, bestMatchBeta);
                else
                    return combinePermutations(bestMatchBeta, bestMatchAlpha, true);
            };


        template<int overallNorm = Eigen::Infinity, int positionalNorm = 2>
        Result compare(ElectronsVector permutee, const ElectronsVector &reference,
                       bool spinSpecificQ = false, bool flipSpinsQ = false) {

            Eigen::PermutationMatrix<Eigen::Dynamic> perm(reference.numberOfEntities());

            if(!spinSpecificQ)
                perm = findPermutation<positionalNorm>(permutee.positionsVector(), reference.positionsVector());
            else
                perm = findSpinSpecificPermutation<positionalNorm>(permutee, reference, flipSpinsQ);

            permutee.permute(perm);

            return {Metrics::positionalNormsVectorNorm<overallNorm, positionalNorm>(
                    permutee.positionsVector(),
                    reference.positionsVector()),
                    std::move(perm)};
        }
    }
}

#endif //INPSIGHTS_BESTMATCHDISTANCE_H
