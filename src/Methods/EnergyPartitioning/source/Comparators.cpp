//
// Created by Michael Heuer on 31.08.18.
//

#include "Comparators.h"
#include <Metrics.h>
#include <Hungarian.h>
#include <Comparators.h>
#include "iomanip"

namespace Comparators {
    ValueEuclideanDistanceComparator::ValueEuclideanDistanceComparator(double distThreshold, double valueThreshold)
    : distThreshold_(distThreshold), valueThreshold_(valueThreshold)
    {}

    bool ValueEuclideanDistanceComparator::operator()(const Reference &lhs, const Reference &rhs) const {

        auto deltaValue = std::abs(lhs.negLogSqrdProbabilityDensity_ - rhs.negLogSqrdProbabilityDensity_);

        if (deltaValue > valueThreshold_) {
            return false;
        }
        else if (deltaValue <= valueThreshold_) {
            return identicalQ(lhs, rhs,distThreshold_);
        } else
            return false;
    }

    bool identicalQ(const Reference &lhs, const Reference &rhs, double distThreshold) {
        assert(lhs.maximum_.typesVector() == rhs.maximum_.typesVector()
               && "The typesvectors must be identical as we assume ordered alpha and beta electrons.");
        assert(lhs.maximum_.positionsVector().numberOfEntities() == rhs.maximum_.positionsVector().numberOfEntities()
               && "The number of positions must be identical.");

        auto lhsCopy = lhs.maximum_.positionsVector();
        auto rhsCopy = rhs.maximum_.positionsVector();

        auto nAlpha = lhs.maximum_.typesVector().countOccurence(Spin::alpha);
        auto nBeta = lhs.maximum_.typesVector().countOccurence(Spin::beta);
        Interval alphaElectrons({0, nAlpha}), betaElectrons({nAlpha, nBeta});

        auto costMatrixAlpha = Metrics::positionalDistances(
                PositionsVector(lhsCopy.slice(alphaElectrons).dataRef()),
                PositionsVector(rhsCopy.slice(alphaElectrons).dataRef()));
        auto costMatrixBeta = Metrics::positionalDistances(
                PositionsVector(lhsCopy.slice(betaElectrons).dataRef()),
                PositionsVector(rhsCopy.slice(betaElectrons).dataRef()));

        auto bestMatchAlpha = Hungarian<double>::findMatching(costMatrixAlpha);
        auto bestMatchBeta = Hungarian<double>::findMatching(costMatrixBeta);

        rhsCopy.slice(alphaElectrons).permute(bestMatchAlpha.inverse()); //TODO INVERSE?
        rhsCopy.slice(betaElectrons).permute(bestMatchBeta.inverse());

        auto distance = Metrics::distance(lhsCopy, rhsCopy);
        return distance <= distThreshold;
    };

    bool globallySimilarQ(const Reference &lhs, const Reference &rhs, double distThreshold) {
        assert(lhs.maximum_.typesVector() == rhs.maximum_.typesVector()
               && "The typesvectors must be identical as we assume ordered alpha and beta electrons.");
        assert(lhs.maximum_.positionsVector().numberOfEntities() == rhs.maximum_.positionsVector().numberOfEntities()
               && "The number of positions must be identical.");

        auto costMatrix = Metrics::positionalDistances(lhs.maximum_.positionsVector(), rhs.maximum_.positionsVector());
        auto bestMatch = Hungarian<double>::findMatching(costMatrix);

        auto rhsCopy = rhs.maximum_;
        rhsCopy.permute(bestMatch.inverse()); //TODO INVERSE?

        auto distance = Metrics::distance(lhs.maximum_.positionsVector(), rhsCopy.positionsVector());
        std::cout << distance << std::endl;

        return distance <= distThreshold;
    };
}