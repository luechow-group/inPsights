

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
        samples[id_].permute(perm);

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


class SimilarReferences {
public:
    explicit SimilarReferences(std::vector<Reference>::iterator representativeReference)
    :
    repRefIt_(representativeReference),
    similarReferences_()
    {}

    void add(std::vector<Reference>::iterator& ref){
        similarReferences_.emplace_back(ref);
    }


    void permuteAll(const Eigen::PermutationMatrix<Eigen::Dynamic>& perm, std::vector<Sample>& samples) {
        (*repRefIt_).permute(perm,samples);
        for (auto& ref : similarReferences_){
            (*ref).permute(perm,samples);
        }
    }

    //TODO REPLACE THIS BY CENTROID LIKE REF
    std::vector<Reference>::iterator repRefIt_; // may change over time, difficult to define for rings/clusters
    std::vector<std::vector<Reference>::iterator> similarReferences_;

};

#endif //AMOLQCPP_REFERENCE_H

