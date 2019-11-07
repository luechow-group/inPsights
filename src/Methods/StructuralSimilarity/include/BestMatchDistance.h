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
