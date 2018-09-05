

//
// Created by Michael Heuer on 28.08.18.
//

#ifndef AMOLQCPP_REFERENCE_H
#define AMOLQCPP_REFERENCE_H

#include <ParticlesVector.h>
#include <Logger.h>

class Reference{
public:
    Reference(double negLogSqrdProbabilityDensity, ElectronsVector maximum = {},  size_t id = 0)
    :
    negLogSqrdProbabilityDensity_(negLogSqrdProbabilityDensity),
    maximum_(std::move(maximum)),
    id_(id),
    associations_({}) // identity perm
    {}

    bool addAssociation(const size_t& id){
        assert(id != id_ && "Self-associations are not allowed");

        if (id != id_) {
            spdlog::get(Logger::name)->warn("Self-associations are not allowed");
            return false;
        }

        auto addedQ = associations_.emplace(id);
        if (!addedQ.second) spdlog::get(Logger::name)->warn("Id is already present.");
        return addedQ.second;
    }

    bool operator <(const Reference& rhs) const {
        return negLogSqrdProbabilityDensity_<rhs.negLogSqrdProbabilityDensity_;
    }

    double negLogSqrdProbabilityDensity_;
    ElectronsVector maximum_;
    size_t id_;
    std::set<size_t> associations_;
};

#endif //AMOLQCPP_REFERENCE_H
