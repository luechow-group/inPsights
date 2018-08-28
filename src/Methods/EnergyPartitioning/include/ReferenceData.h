#include <utility>

//
// Created by Michael Heuer on 28.08.18.
//

#ifndef AMOLQCPP_REFERENCEDATA_H
#define AMOLQCPP_REFERENCEDATA_H

#include <map>
#include <ParticlesVector.h>
#include "SampleData.h"



class Reference{
public:

    Reference(ElectronsVector maximum, double negLogSqrdProbabilityDensity)
    : maximum_(std::move(maximum)), negLogSqrdProbabilityDensity_(negLogSqrdProbabilityDensity)
    {}

    bool operator<(const Reference& other) const {
        if (std::abs(negLogSqrdProbabilityDensity_-other.negLogSqrdProbabilityDensity_) > valueThreshold)
            return true;
        else if (std::abs(negLogSqrdProbabilityDensity_-other.negLogSqrdProbabilityDensity_) <= valueThreshold)
            return true; //TODO hungarian + euclidean distance check
        else
            return false;
    }


    ElectronsVector maximum_;
    double negLogSqrdProbabilityDensity_;

    const double distThreshold = 0.01;
    const double valueThreshold = 1e-4;
};

using RefSamplePair = std::pair<Reference,Sample>;
class ReferenceSampleMapping {
public:
    std::multimap<Reference,Sample> map;
};




#endif //AMOLQCPP_REFERENCEDATA_H
