//
// Created by heuer on 19.10.18.
//

#ifndef AMOLQCPP_SIMILARREFERENCES_H
#define AMOLQCPP_SIMILARREFERENCES_H

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

#endif //AMOLQCPP_SIMILARREFERENCES_H
