

//
// Created by Michael Heuer on 28.08.18.
//

#ifndef AMOLQCPP_REFERENCE_H
#define AMOLQCPP_REFERENCE_H

#include <ParticlesVector.h>
#include <Sample.h>
#include <Statistics.h>

class Reference{
public:
    Reference(double negLogSqrdProbabilityDensity, ElectronsVector maximum = {}, size_t id = 0)
    :
    negLogSqrdProbabilityDensity_(negLogSqrdProbabilityDensity),
    maximum_(std::move(maximum)),
    sampleIds_({id})
    {}

    size_t ownId() const {
        return sampleIds_[0];
    }

    //TODO is this the task of a container?
    void mergeReference(std::vector<Reference>::iterator &it) {
        assert((*it).ownId() != ownId() && "Self-associations are not allowed");

        sampleIds_.insert(
                sampleIds_.end(),
                make_move_iterator((*it).sampleIds_.begin()),
                make_move_iterator((*it).sampleIds_.end()));
    }

    void permute(const Eigen::PermutationMatrix<Eigen::Dynamic>& perm, std::vector<Sample>& samples){
        maximum_.permute(perm);

        for(const auto & i : sampleIds_){
            samples[i].permute(perm);
        }
    }

    bool operator<(const Reference& rhs) const {
        return negLogSqrdProbabilityDensity_<rhs.negLogSqrdProbabilityDensity_;
    }

    unsigned long count() const {
        return sampleIds_.size();
    }

    double value() const {
        return negLogSqrdProbabilityDensity_;
    }

    ElectronsVector maximum() const {
        return maximum_;
    }

    std::vector<size_t> sampleIds() const {
        return sampleIds_;
    }


private:
    double negLogSqrdProbabilityDensity_;
    ElectronsVector maximum_;
    std::vector<size_t> sampleIds_; // contain samples instead?
};



#endif //AMOLQCPP_REFERENCE_H

