//
// Created by Michael Heuer on 31.08.18.
//

#include "Comparators.h"
#include <Metrics.h>
#include <Hungarian.h>

bool ValueEuclideanDistanceComparator::operator() (const Reference& lhs, const Reference& rhs) const {

    auto deltaValue = std::abs(lhs.negLogSqrdProbabilityDensity_-rhs.negLogSqrdProbabilityDensity_);

    std::cout << deltaValue << std::endl;

    if (deltaValue > valueThreshold)
        return false;
    else if (deltaValue <= valueThreshold) {
        auto costMatrix = Metrics::positionalDistances(lhs.maximum_.positionsVector(),rhs.maximum_.positionsVector());
        auto permMatrix = Hungarian<double>::findMatching(costMatrix);

        auto copy = lhs.maximum_;
        copy.permute(permMatrix);

        auto distance = Metrics::distance(copy.positionsVector(),rhs.maximum_.positionsVector());
        std::cout << distance<< std::endl;

        return distance <= distThreshold;
    }
    else
        return false;
};