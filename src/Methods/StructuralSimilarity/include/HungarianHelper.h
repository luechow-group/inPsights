//
// Created by Michael Heuer on 06.09.18.
//

#ifndef INPSIGHTS_HUNGARIANHELPER_H
#define INPSIGHTS_HUNGARIANHELPER_H


#include "Hungarian.h"
#include <Metrics.h>
#include <ParticlesVector.h>
#include <Interval.h>
#include <StructuralSimilarity.h>

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


    template<int overallNorm = Eigen::Infinity, int positionalNorm = 2>
    std::pair<double,Eigen::PermutationMatrix<Eigen::Dynamic>> bestMatchSimilarity(
            const MolecularSpectrum &permutee,
            const MolecularSpectrum &reference
            ) {
        assert(ParticleKit::isSubsetQ(permutee.molecule_) && "The permutee must be a subset of the particle kit.");
        assert(ParticleKit::isSubsetQ(reference.molecule_) && "The reference must be a subset of the particle kit.");

        auto costMatrix = StructuralSimilarity::correlationMatrix(permutee, reference);
        Eigen::PermutationMatrix<Eigen::Dynamic> bestMatch = Hungarian<double>::findMatching(costMatrix, Matchtype::MAX);

        // best-match permute columns and sum diagonal elements
        double simMetric = (costMatrix * bestMatch).diagonal().sum();
        return {simMetric, std::move(bestMatch)};
    }
}

#endif //INPSIGHTS_HUNGARIANHELPER_H
