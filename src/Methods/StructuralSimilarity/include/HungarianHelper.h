//
// Created by Michael Heuer on 06.09.18.
//

#ifndef AMOLQCPP_HUNGARIANHELPER_H
#define AMOLQCPP_HUNGARIANHELPER_H


#include "Hungarian.h"
#include <Metrics.h>
#include <ParticlesVector.h>

namespace HungarianHelper{

    Eigen::PermutationMatrix<Eigen::Dynamic> combinePermutations(
            const Eigen::PermutationMatrix<Eigen::Dynamic>& p1,
            const Eigen::PermutationMatrix<Eigen::Dynamic>& p2) {

        long n1 = p1.size(), n2 = p2.size();
        Eigen::VectorXi combined(n1+n2);

        combined.segment(0,n1) = p1.indices().base();
        combined.segment(n1,n2) = (p2.indices().base().array()+n1);

        return Eigen::PermutationMatrix<Eigen::Dynamic>(combined);
    };

    Eigen::PermutationMatrix<Eigen::Dynamic> spinSpecificHungarian(
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
            // flip slice intervals
            rhsAlpha = lhsBeta;
            rhsBeta = lhsAlpha;
        }

        auto costMatrixAlpha = Metrics::positionalDistances(
                PositionsVector(lhsCopy.slice(lhsAlpha).dataRef()),
                PositionsVector(rhsCopy.slice(rhsAlpha).dataRef()));
        auto costMatrixBeta = Metrics::positionalDistances(
                PositionsVector(lhsCopy.slice(lhsBeta).dataRef()),
                PositionsVector(rhsCopy.slice(rhsBeta).dataRef()));

        auto bestMatchAlpha = Hungarian<double>::findMatching(costMatrixAlpha);
        auto bestMatchBeta = Hungarian<double>::findMatching(costMatrixBeta);

        return combinePermutations(bestMatchAlpha,bestMatchBeta);
    };

    double mostDeviatingParticleDistance(const PositionsVector& ref, PositionsVector other,
            const Eigen::PermutationMatrix<Eigen::Dynamic>& perm) {
        assert(ref.numberOfEntities() == other.numberOfEntities());

        other.permute(perm);
        return Metrics::positionDistancesVector(ref,other).lpNorm<Eigen::Infinity>();
    }
};

#endif //AMOLQCPP_HUNGARIANHELPER_H
