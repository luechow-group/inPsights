//
// Created by Michael Heuer on 06.09.18.
//

#ifndef INPSIGHTS_BESTMATCH_H
#define INPSIGHTS_BESTMATCH_H

#include "Hungarian.h"
#include <Metrics.h>
#include <ParticlesVector.h>
#include <Interval.h>
#include <LocalSimilarity.h>
#include <MolecularSpectrum.h>
#include <ParticleKit.h>

namespace BestMatch {
    Eigen::PermutationMatrix<Eigen::Dynamic> combinePermutations(
            const Eigen::PermutationMatrix<Eigen::Dynamic> &p1,
            const Eigen::PermutationMatrix<Eigen::Dynamic> &p2, bool flipSpinsQ = false);

    struct Result {
        const double metric;
        const Eigen::PermutationMatrix<Eigen::Dynamic> permutation;
    };

    namespace Distance {
        template<int positionalNorm = 2>
        Eigen::PermutationMatrix<Eigen::Dynamic> findPermutation(
                const PositionsVector &permutee,
                const PositionsVector &reference) {
            assert(permutee.numberOfEntities() == reference.numberOfEntities());

            auto costMatrix = Metrics::positionalDistances<positionalNorm>(permutee, reference);

            return Hungarian<double>::findMatching(costMatrix);
        }


        template<int overallNorm = Eigen::Infinity, int positionalNorm = 2>
        Result compare(ElectronsVector permutee, const ElectronsVector &reference) {
            auto perm = findPermutation<positionalNorm>(permutee.positionsVector(), reference.positionsVector());

            permutee.permute(perm);

            return {Metrics::positionalNormsVectorNorm<overallNorm, positionalNorm>(
                    permutee.positionsVector(),
                    reference.positionsVector()),
                    std::move(perm)};
        }


        namespace SpinSpecific {
            template<int positionalNorm = 2>//TODO !!!!!!!!!!! Check if rhs is reference and lhs is permutee
            Eigen::PermutationMatrix<Eigen::Dynamic> findPermutation(ElectronsVector lhs, ElectronsVector rhs, bool flipSpinsQ = false) {
                assert(lhs.typesVector() == rhs.typesVector()
                       && "The typesvectors must be identical as we assume ordered alpha and beta electrons.");
                assert(lhs.positionsVector().numberOfEntities() == rhs.positionsVector().numberOfEntities()
                       && "The number of positions must be identical.");

                auto nAlpha = lhs.typesVector().countOccurence(Spin::alpha);
                auto nBeta = lhs.typesVector().countOccurence(Spin::beta);

                Interval
                        lhsAlpha{0, nAlpha}, rhsAlpha{},
                        lhsBeta{nAlpha, nBeta}, rhsBeta{};

                if (!flipSpinsQ) {
                    rhsAlpha = lhsAlpha;
                    rhsBeta = lhsBeta;
                } else {
                    assert(nAlpha == nBeta
                    && "The number of alpha and beta electrons must match to allow for spin flips.");

                    // flip slice intervals
                    rhsAlpha = lhsBeta;
                    rhsBeta = lhsAlpha;
                };

                auto costMatrixAlpha = Metrics::positionalDistances<positionalNorm>(
                        PositionsVector(lhs.positionsVector().asEigenVector().segment(lhsAlpha.start() * 3,
                                        lhsAlpha.numberOfEntities() * 3)),
                        PositionsVector(
                                rhs.positionsVector().asEigenVector().segment(rhsAlpha.start() * 3,
                                                                rhsAlpha.numberOfEntities() * 3)));

                auto bestMatchAlpha = Hungarian<double>::findMatching(costMatrixAlpha);

                auto costMatrixBeta = Metrics::positionalDistances<positionalNorm>(
                        PositionsVector(lhs.positionsVector().asEigenVector().segment(lhsBeta.start() * 3,
                                        lhsBeta.numberOfEntities() * 3)),
                        PositionsVector(rhs.positionsVector().asEigenVector().segment(rhsBeta.start() * 3,
                                        rhsBeta.numberOfEntities() * 3)));

                auto bestMatchBeta = Hungarian<double>::findMatching(costMatrixBeta);

                if (!flipSpinsQ)
                    return combinePermutations(bestMatchAlpha, bestMatchBeta);
                else
                    return combinePermutations(bestMatchBeta, bestMatchAlpha, true);
            };

            template<int overallNorm = Eigen::Infinity, int positionalNorm = 2>
            Result compare(ElectronsVector permutee, const ElectronsVector &reference, bool flipSpinsQ = false) {

                auto perm = findPermutation<positionalNorm>(permutee, reference, flipSpinsQ);
                permutee.permute(perm);

                return {Metrics::positionalNormsVectorNorm<overallNorm, positionalNorm>(
                        permutee.positionsVector(), reference.positionsVector()),
                        std::move(perm)
                };
            }
        }
    }

    namespace Similarity {
        Result compare(
                const MolecularSpectrum &permutee,
                const MolecularSpectrum &reference);
    }
};

#endif //INPSIGHTS_BESTMATCH_H
