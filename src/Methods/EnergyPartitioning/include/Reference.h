

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


private:
    double negLogSqrdProbabilityDensity_;
    ElectronsVector maximum_;
    std::vector<size_t> sampleIds_; // contain samples instead?
};


class SimilarReferences {
public:
    explicit SimilarReferences(std::vector<Reference>::iterator representativeReference)
    :
    similarReferences_({representativeReference}),
    valueStats_()
    {
        Eigen::Matrix<double,1,1> v(1);
        v << (*representativeReference).value();
        valueStats_.add(v);
    }

    void add(std::vector<Reference>::iterator& ref){
        similarReferences_.emplace_back(ref);
        Eigen::Matrix<double,1,1> v(1);
        v << (*ref).value();
        valueStats_.add(v);
    }

    const Reference& representativeReference() const {
        return *similarReferences_[0];
    }

    Reference& representativeReference(){
        return *similarReferences_[0];
    }

    std::vector<Reference>::iterator representativeReferenceIterator(){
        return similarReferences_[0];
    }

    void permuteAll(const Eigen::PermutationMatrix<Eigen::Dynamic>& perm, std::vector<Sample>& samples) {
        for (auto& ref : similarReferences_){
            (*ref).permute(perm,samples);
        }
    }

    std::vector<std::vector<Reference>::iterator>& similarReferencesIterators(){
        return similarReferences_;
    }

    const std::vector<std::vector<Reference>::iterator>& similarReferencesIterators() const {
        return similarReferences_;
    }

    const Statistics::RunningStatistics<Eigen::Matrix<double,1,1>>& valueStats() const {
        return valueStats_;
    }

    bool operator<(const SimilarReferences& rhs) const {
        //std::cout << valueStats_.cwiseMin()[0] << ", " << rhs.valueStats().cwiseMin()[0] << std::endl;
        return representativeReference().value() < rhs.representativeReference().value();
        //return valueStats_.cwiseMin()[0] < rhs.valueStats().cwiseMin()[0];
    }

    //TODO REPLACE THIS BY CENTROID LIKE REF
private:
    std::vector<std::vector<Reference>::iterator> similarReferences_;
    Statistics::RunningStatistics<Eigen::Matrix<double,1,1>> valueStats_;
};

#endif //AMOLQCPP_REFERENCE_H

