//
// Created by heuer on 03.04.19.
//

#ifndef INPSIGHTS_BESTMATCHDISTANCE_H
#define INPSIGHTS_BESTMATCHDISTANCE_H

#include "Hungarian.h"
#include <BestMatch.h>
#include "BestMatchDistance.h"
#include <Metrics.h>

namespace BestMatch {
    namespace Distance {
        template<int positionalNorm = 2> Eigen::PermutationMatrix<Eigen::Dynamic> findPermutation(
                const PositionsVector &permutee, const PositionsVector &reference) {
            assert(permutee.numberOfEntities() == reference.numberOfEntities());

            auto costMatrix = Metrics::positionalDistances<positionalNorm>(permutee, reference);

            return Hungarian<double>::findMatching(costMatrix);
        }


        template<typename Type, int positionalNorm = 2>
        Eigen::PermutationMatrix<Eigen::Dynamic>
        findTypeSpecificPermutation(const ParticlesVector<Type>& permutee, const ParticlesVector<Type>& reference) {
            assert(permutee.numberOfEntities() == reference.numberOfEntities()
            && "The number of particles must be identical.");

            auto countedTypes = permutee.typesVector().countTypes();
            
            assert(countedTypes == reference.typesVector().countTypes() 
            && "The type vectors must have the same type counts.");

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

        /* Type unspecific comparison. The Type attribute is only necessary to allow the method to work for all kinds
         * of ParticleVectors.
         */
        template<int overallNorm = Eigen::Infinity, int positionalNorm = 2>
        Result compare(PositionsVector permutee, const PositionsVector &reference) {

            auto perm = findPermutation<positionalNorm>(permutee, reference);

            permutee.permute(perm);

            return {Metrics::positionalNormsVectorNorm<overallNorm, positionalNorm>(permutee,reference),
                    std::move(perm)};
        }

        // Type specific comparison
        template<typename Type, int overallNorm = Eigen::Infinity, int positionalNorm = 2>
        Result compare(ParticlesVector<Type> permutee, const ParticlesVector<Type> &reference) {

            auto perm = findTypeSpecificPermutation<Type,positionalNorm>(permutee, reference);

            permutee.permute(perm);

            return {Metrics::positionalNormsVectorNorm<overallNorm, positionalNorm>(
                    permutee.positionsVector(),
                    reference.positionsVector()),
                    std::move(perm)};
        }
    }
}

#endif //INPSIGHTS_BESTMATCHDISTANCE_H
