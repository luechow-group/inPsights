

//
// Created by Michael Heuer on 28.08.18.
//

#ifndef AMOLQCPP_REFERENCEDATA_H
#define AMOLQCPP_REFERENCEDATA_H

#include <ParticlesVector.h>
#include <utility>

class Reference{
public:

    Reference()= default;

    Reference(ElectronsVector maximum, double negLogSqrdProbabilityDensity)
    : maximum_(std::move(maximum)), negLogSqrdProbabilityDensity_(negLogSqrdProbabilityDensity)
    {}

    ElectronsVector maximum_;
    double negLogSqrdProbabilityDensity_;

private:

};






#endif //AMOLQCPP_REFERENCEDATA_H
