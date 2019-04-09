//
// Created by heuer on 19.10.18.
//

#ifndef INPSIGHTS_SIMILARREFERENCES_H
#define INPSIGHTS_SIMILARREFERENCES_H

#include "Reference.h"

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
        for (auto ref : similarReferences_)
            ref.base()->permute(perm,samples);
    }

    std::vector<std::vector<Reference>::iterator>& similarReferencesIterators(){
        return similarReferences_;
    }

    const std::vector<std::vector<Reference>::iterator>& similarReferencesIterators() const {
        return similarReferences_;
    }

    const SingleValueStatistics& valueStats() const {
        return valueStats_;
    }

    bool operator<(const SimilarReferences& rhs) const {
        return representativeReference().value() < rhs.representativeReference().value();
    }

    void sort(){
        std::sort(similarReferences_.begin(), similarReferences_.end());
    }

private:
    std::vector<std::vector<Reference>::iterator> similarReferences_;
    SingleValueStatistics valueStats_;
};

#endif //INPSIGHTS_SIMILARREFERENCES_H
