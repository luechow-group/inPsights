//
// Created by Michael Heuer on 06.09.18.
//

#ifndef INPSIGHTS_HUNGARIANHELPER_H
#define INPSIGHTS_HUNGARIANHELPER_H


#include "Hungarian.h"
#include <Metrics.h>
#include <ParticlesVector.h>
#include <Interval.h>
#include <LocalSimilarity.h>
#include <MolecularSpectrum.h>
#include <ParticleKit.h>

namespace HungarianHelper{

    Eigen::PermutationMatrix<Eigen::Dynamic> combinePermutations(
            const Eigen::PermutationMatrix<Eigen::Dynamic>& p1,
            const Eigen::PermutationMatrix<Eigen::Dynamic>& p2, bool flipSpinsQ = false);



    template <int positionalNorm = 2>
    Eigen::PermutationMatrix<Eigen::Dynamic> spinSpecificBestMatch(
            const ElectronsVector &lhs,
            const ElectronsVector &rhs,
            bool flipSpinsQ = false) {
        assert(lhs.typesVector() == rhs.typesVector()
               && "The typesvectors must be identical as we assume ordered alpha and beta electrons.");
        assert(lhs.positionsVector().numberOfEntities() == rhs.positionsVector().numberOfEntities()
               && "The number of positions must be identical.");

        auto lhsCopy = lhs.positionsVector();//TODO eliminate redundant copy by const_cast .slice
        //https://stackoverflow.com/questions/5008541/how-to-call-a-non-const-function-within-a-const-function-c
        auto rhsCopy = rhs.positionsVector();

        auto nAlpha = lhs.typesVector().countOccurence(Spin::alpha);
        auto nBeta = lhs.typesVector().countOccurence(Spin::beta);

        Interval
        lhsAlpha{0, nAlpha}, rhsAlpha{},
        lhsBeta{nAlpha, nBeta}, rhsBeta{};

        if(!flipSpinsQ){
            rhsAlpha = lhsAlpha;
            rhsBeta = lhsBeta;
        } else {
            assert(nAlpha == nBeta && "The number of alpha and beta electrons must match to allow for spin flips.");
            // flip slice intervals
            rhsAlpha = lhsBeta;
            rhsBeta = lhsAlpha;
        };

        auto costMatrixAlpha = Metrics::positionalDistances<positionalNorm>(
                PositionsVector(lhsCopy.asEigenVector().segment(lhsAlpha.start()*3,lhsAlpha.numberOfEntities()*3)),
                PositionsVector(rhsCopy.asEigenVector().segment(rhsAlpha.start()*3,rhsAlpha.numberOfEntities()*3)));
        auto bestMatchAlpha = Hungarian<double>::findMatching(costMatrixAlpha);

        auto costMatrixBeta = Metrics::positionalDistances<positionalNorm>(
                PositionsVector(lhsCopy.asEigenVector().segment(lhsBeta.start()*3,lhsBeta.numberOfEntities()*3)),
                PositionsVector(rhsCopy.asEigenVector().segment(rhsBeta.start()*3,rhsBeta.numberOfEntities()*3)));
        auto bestMatchBeta = Hungarian<double>::findMatching(costMatrixBeta);

        if(!flipSpinsQ)
            return combinePermutations(bestMatchAlpha, bestMatchBeta);
        else
            return combinePermutations(bestMatchBeta, bestMatchAlpha, true);
    };
};

namespace Metrics{
    template<int overallNorm = Eigen::Infinity, int positionalNorm = 2>
    double bestMatchNorm(
            const PositionsVector &permutee,
            const Eigen::PermutationMatrix<Eigen::Dynamic> &perm,
            const PositionsVector &reference) {
        assert(permutee.numberOfEntities() == reference.numberOfEntities());

        auto copy = permutee;
        copy.permute(perm);

        Eigen::VectorXd res = Metrics::positionalNormsVector<positionalNorm>(copy, reference);
        return res.lpNorm<overallNorm>();
    }

    template<int overallNorm = Eigen::Infinity, int positionalNorm = 2>
    std::pair<double,Eigen::PermutationMatrix<Eigen::Dynamic>> bestMatch(
            const PositionsVector &permutee,
            const PositionsVector &reference) {
        assert(permutee.numberOfEntities() == reference.numberOfEntities());

        auto costMatrix = Metrics::positionalDistances<positionalNorm>(permutee,reference);
        Eigen::PermutationMatrix<Eigen::Dynamic> bestMatch = Hungarian<double>::findMatching(costMatrix);

        return {bestMatchNorm<overallNorm, positionalNorm>(permutee, bestMatch, reference), std::move(bestMatch)};
    }

    template<int overallNorm = Eigen::Infinity, int positionalNorm = 2>
    double bestMatchNorm(const PositionsVector &permutee, const PositionsVector &reference) {
        return bestMatch<overallNorm,positionalNorm>(permutee,reference).first;
    }

    template<int overallNorm = Eigen::Infinity, int positionalNorm = 2>
    std::pair<double,Eigen::PermutationMatrix<Eigen::Dynamic>> bestMatch(
            const ElectronsVector &permutee,
            const ElectronsVector &reference){
        return bestMatch<overallNorm, positionalNorm>(permutee.positionsVector(), reference.positionsVector());
    }

    template<int overallNorm = Eigen::Infinity, int positionalNorm = 2>
    double bestMatchNorm(const ElectronsVector &permutee, const ElectronsVector &reference){
        return bestMatch<overallNorm,positionalNorm>(permutee,reference).first;
    }

    template<int overallNorm = Eigen::Infinity, int positionalNorm = 2>
    std::pair<double,Eigen::PermutationMatrix<Eigen::Dynamic>> spinSpecificBestMatch(
            const ElectronsVector &permutee,
            const ElectronsVector &reference,
            bool flipSpinsQ = false){

        auto spinSpecificBestMatch = HungarianHelper::spinSpecificBestMatch<positionalNorm>(permutee,reference,flipSpinsQ);

        return {bestMatchNorm<overallNorm, positionalNorm>(
                permutee.positionsVector(), spinSpecificBestMatch,
                reference.positionsVector()),
                std::move(spinSpecificBestMatch)};
    }

    template<int overallNorm = Eigen::Infinity, int positionalNorm = 2>
    double spinSpecificBestMatchNorm(
            const ElectronsVector &permutee, const ElectronsVector &reference, bool flipSpinsQ = false){
        return spinSpecificBestMatch<overallNorm,positionalNorm>(permutee,reference,flipSpinsQ).first;
    }


    std::pair<double,Eigen::PermutationMatrix<Eigen::Dynamic>> bestMatchSimilarity(
            const MolecularSpectrum &permutee,
            const MolecularSpectrum &reference) {
        assert(ParticleKit::isSubsetQ(permutee.molecule_) && "The permutee must be a subset of the particle kit.");
        assert(ParticleKit::isSubsetQ(reference.molecule_) && "The reference must be a subset of the particle kit.");

        // TODO assert that identical number of electrons and same atom geometry? Is this constraint needed? What happens with rows/cols of zero?

        assert(
                permutee.molecule_.electrons().typesVector().countOccurence(Spin::alpha) ==
                reference.molecule_.electrons().typesVector().countOccurence(Spin::alpha)
                && "The number of alpha electrons has to match.");

        assert(
                permutee.molecule_.electrons().typesVector().countOccurence(Spin::beta) ==
                reference.molecule_.electrons().typesVector().countOccurence(Spin::beta)
                && "The number of beta electrons has to match.");

        auto nAlpha = reference.molecule_.electrons().typesVector().countOccurence(Spin::alpha);
        auto nBeta = reference.molecule_.electrons().typesVector().countOccurence(Spin::beta);

        auto N = nAlpha+nBeta;
        Eigen::MatrixXd environmentalSimilarities(N,N);

        // TODO consider identical spin flip
        TypeSpecificNeighborhoodsAtOneCenter expA, expB;
        for (unsigned i = 0; i < nAlpha; ++i) {
            EnumeratedType<int> enumeratedType_i(Spins::spinToInt(Spin::alpha), i);
            expA = permutee.molecularCenters_.find(enumeratedType_i)->second;

            for (unsigned j = 0; j < nAlpha; ++j) {
                EnumeratedType<int> enumeratedType_j(Spins::spinToInt(Spin::alpha), j);
                expB = reference.molecularCenters_.find(enumeratedType_j)->second;
                environmentalSimilarities(i, j) = LocalSimilarity::kernel(expA, expB,SOAPExpansion::settings.zeta());
            }
            for (unsigned j = 0; j < nBeta; ++j) {
                EnumeratedType<int> enumeratedType_j(Spins::spinToInt(Spin::beta), j);
                expB = reference.molecularCenters_.find(enumeratedType_j)->second;
                environmentalSimilarities(i, nAlpha+j) = LocalSimilarity::kernel(expA, expB,SOAPExpansion::settings.zeta());
            }
        }
        for (unsigned i = 0; i < nBeta; ++i) {
            EnumeratedType<int> enumeratedType_i(Spins::spinToInt(Spin::beta), i);
            expA = permutee.molecularCenters_.find(enumeratedType_i)->second;

            for (unsigned j = 0; j < nAlpha; ++j) {
                EnumeratedType<int> enumeratedType_j(Spins::spinToInt(Spin::alpha), j);
                expB = reference.molecularCenters_.find(enumeratedType_j)->second;
                environmentalSimilarities(nAlpha+i, j) = LocalSimilarity::kernel(expA, expB,SOAPExpansion::settings.zeta());
            }
            for (unsigned j = 0; j < nBeta; ++j) {
                EnumeratedType<int> enumeratedType_j(Spins::spinToInt(Spin::beta), j);
                expB = reference.molecularCenters_.find(enumeratedType_j)->second;
                environmentalSimilarities(nAlpha+i, nAlpha+j) = LocalSimilarity::kernel(expA, expB,SOAPExpansion::settings.zeta());
            }
        }

        Eigen::PermutationMatrix<Eigen::Dynamic> bestMatch = Hungarian<double>::findMatching(environmentalSimilarities, Matchtype::MAX);

        // best-match permute columns and sum diagonal elements
        double simMetric = (environmentalSimilarities * bestMatch).diagonal().sum() / N;

        //restore the original order before the particle kit permutations
        auto permuteeToKit = ParticleKit::toKitPermutation(permutee.molecule_.electrons());
        auto referenceFromKit= ParticleKit::fromKitPermutation(reference.molecule_.electrons());

        return {simMetric, referenceFromKit * bestMatch * permuteeToKit};
    }
}

#endif //INPSIGHTS_HUNGARIANHELPER_H
