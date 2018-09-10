

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
    associatedSampleIds_({})
    {}

    void addAssociation(const size_t& id) {
        assert(id != id_ && "Self-associations are not allowed");
        associatedSampleIds_.emplace_back(id);
    }

    bool operator<(const Reference& rhs) const {
        return negLogSqrdProbabilityDensity_<rhs.negLogSqrdProbabilityDensity_;
    }

    double negLogSqrdProbabilityDensity_;
    ElectronsVector maximum_;
    size_t id_;
    std::vector<size_t> associatedSampleIds_; //associated Samples id/identical Maxima
};


class SimilarReference {
public:
    SimilarReference(std::vector<Reference>::iterator ref, const Eigen::PermutationMatrix<Eigen::Dynamic> &perm)
    : it_(ref), perm_(perm) {}

    bool operator <(const SimilarReference& rhs) const {
        return (*it_).negLogSqrdProbabilityDensity_< (*rhs.it_).negLogSqrdProbabilityDensity_;
    }

    std::vector<Reference>::iterator it_; // association to the list of globally identical maxima
    Eigen::PermutationMatrix<Eigen::Dynamic> perm_; // perm to turn
};


class SimilarReferencesCollection {
public:
    explicit SimilarReferencesCollection(std::vector<Reference>::iterator representativeReference)
    : representativeReferenceIterator(representativeReference) {}
    //TODO REPLACE THIS BY CENTROID LIKE REF

    std::vector<Reference>::iterator representativeReferenceIterator; // may change over time, difficult to define for rings/clusters
    std::vector<SimilarReference> similarReferences_;

};

#endif //AMOLQCPP_REFERENCE_H

