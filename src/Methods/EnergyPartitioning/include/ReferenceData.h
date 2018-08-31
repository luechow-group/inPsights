#include <utility>

//
// Created by Michael Heuer on 28.08.18.
//

#ifndef AMOLQCPP_REFERENCEDATA_H
#define AMOLQCPP_REFERENCEDATA_H

#include <map>
#include <ParticlesVector.h>
#include "SampleData.h"
#include "Metrics.h"
#include <Hungarian.h>

class Reference{
public:

    Reference(ElectronsVector maximum, double negLogSqrdProbabilityDensity)
    : maximum_(std::move(maximum)), negLogSqrdProbabilityDensity_(negLogSqrdProbabilityDensity)
    {}

    bool operator<(const Reference& other) const {
        if (std::abs(negLogSqrdProbabilityDensity_-other.negLogSqrdProbabilityDensity_) > valueThreshold)
            return true;
        else if (std::abs(negLogSqrdProbabilityDensity_-other.negLogSqrdProbabilityDensity_) <= valueThreshold) {

            auto costMatrix = Metrics::positionalDistances(maximum_.positionsVector(),other.maximum_.positionsVector());
            auto permMatrix = Hungarian<double>::findMatching(costMatrix);

            auto copy = maximum_;
            copy.permute(permMatrix);

            auto distance = Metrics::distance(copy.positionsVector(),other.maximum_.positionsVector());

            return distance <= distThreshold;
        }
        else
            return false;
    }


    ElectronsVector maximum_;
    double negLogSqrdProbabilityDensity_;

private:
    const double distThreshold = 0.01;
    const double valueThreshold = 1e-4;
};

using RefSamplePair = std::pair<Reference,Sample>;

class ReferenceSampleMapping {
public:
    std::multimap<Reference,Sample> map;
};




#endif //AMOLQCPP_REFERENCEDATA_H
