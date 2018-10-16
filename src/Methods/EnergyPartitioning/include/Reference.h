

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
    ids_({id})
    {}

    size_t ownId() const {
        return ids_[0];
    }

    void addSampleIds(std::vector<Reference>::iterator &it) {
        assert((*it).ownId() != ownId() && "Self-associations are not allowed");

        ids_.insert(
                ids_.end(),
                make_move_iterator((*it).ids_.begin()),
                make_move_iterator((*it).ids_.end()));
    }

    void permute(const Eigen::PermutationMatrix<Eigen::Dynamic>& perm, std::vector<Sample>& samples){
        maximum_.permute(perm);

        for(const auto & i : ids_){
            samples[i].permute(perm);
        }
    }

    bool operator<(const Reference& rhs) const {
        return negLogSqrdProbabilityDensity_<rhs.negLogSqrdProbabilityDensity_;
    }

    unsigned long count() const {
        return ids_.size();
    }

    double negLogSqrdProbabilityDensity_;
    ElectronsVector maximum_;
    std::vector<size_t> ids_;
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

