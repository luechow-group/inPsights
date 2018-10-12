

//
// Created by Michael Heuer on 28.08.18.
//

#ifndef AMOLQCPP_REFERENCE_H
#define AMOLQCPP_REFERENCE_H

#include <ParticlesVector.h>
#include <Sample.h>

class Reference{
public:
    Reference(double negLogSqrdProbabilityDensity, ElectronsVector maximum = {},  size_t id = 0)
    :
    negLogSqrdProbabilityDensity_(negLogSqrdProbabilityDensity),
    maximum_(std::move(maximum)),
    id_(id),
    associatedSampleIds_({})
    {}

    void addAssociations(std::vector<Reference>::iterator& it) {
        assert((*it).id_ != id_ && "Self-associations are not allowed");

        associatedSampleIds_.emplace_back((*it).id_);

        associatedSampleIds_.insert(
                associatedSampleIds_.end(),
                make_move_iterator((*it).associatedSampleIds_.begin()),
                make_move_iterator((*it).associatedSampleIds_.end()));
    }

    void permute(const Eigen::PermutationMatrix<Eigen::Dynamic>& perm, std::vector<Sample>& samples){
        maximum_.permute(perm);
        for(const auto & i : associatedSampleIds_){
            samples[i].permute(perm);
        }
    }

    bool operator<(const Reference& rhs) const {
        return negLogSqrdProbabilityDensity_<rhs.negLogSqrdProbabilityDensity_;
    }

    unsigned long count() const {
        return 1+associatedSampleIds_.size();
    }

    double negLogSqrdProbabilityDensity_;
    ElectronsVector maximum_;
    size_t id_;
    std::vector<size_t> associatedSampleIds_; //associated Samples id/identical Maxima
};


class SimilarReference { //TODO delete class.
public:
    SimilarReference(std::vector<Reference>::iterator ref)
    : it_(ref) {}

    bool operator <(const SimilarReference& rhs) const {
        return (*it_).negLogSqrdProbabilityDensity_< (*rhs.it_).negLogSqrdProbabilityDensity_;
    }

    std::vector<Reference>::iterator it_; // association to the list of globally identical maxima
};


class SimilarReferences {
public:
    explicit SimilarReferences(std::vector<Reference>::iterator representativeReference)
    :
    repRefIt_(representativeReference),
    similarReferences_()
    {}

    //TODO REPLACE THIS BY CENTROID LIKE REF
    std::vector<Reference>::iterator repRefIt_; // may change over time, difficult to define for rings/clusters
    std::vector<SimilarReference> similarReferences_;

};

#endif //AMOLQCPP_REFERENCE_H

