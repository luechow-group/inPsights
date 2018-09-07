

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

    bool addAssociation(const size_t& id) const {//TODO BE CAREFUL WITH CONST
        assert(id != id_ && "Self-associations are not allowed");

        if (id == id_) {
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
    mutable std::set<size_t> associations_;
};


class SimilarReference {
public:
    SimilarReference(std::unique_ptr<Reference> ref, Eigen::PermutationMatrix<Eigen::Dynamic> perm)
            :
            ref_(std::move(ref)), perm_(std::move(perm))
    {}

    std::unique_ptr<Reference> ref_;
    Eigen::PermutationMatrix<Eigen::Dynamic> perm_; // perm to turn
};


class SimilarReferencesCollection {
public:
    /*TODO CONTINUE HERE!!!!!!!!!!!!*/
    SimilarReferencesCollection()//std::unique_ptr<Reference> representativeReference)
    :
    representativeReference_(std::move(representativeReference))
    {}

    std::unique_ptr<Reference> representativeReference_; // may change over time, difficult to define for rings/clusters
    std::set<SimilarReference> similarReferences_;

};



#endif //AMOLQCPP_REFERENCE_H
