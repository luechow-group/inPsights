//
// Created by Michael Heuer on 31.08.18.
//

#ifndef INPSIGHTS_DATAVECTOR_H
#define INPSIGHTS_DATAVECTOR_H

#include "AbstractVector.h"
#include <Eigen/Core>

template<typename Scalar>
class DataVector : public AbstractVector {

protected:
    using EigenVecType = Eigen::Matrix<Scalar, Eigen::Dynamic,1>;
    using RefEigenVecType = Eigen::Ref<EigenVecType>;

public:
    DataVector(long numberOfEntities = 0, long entityLength = 1)
    :
    AbstractVector(numberOfEntities,entityLength),
    data_(EigenVecType::Constant(numberOfEntities, 0))
    {}

    DataVector(const DataVector &rhs)
    :
    AbstractVector(rhs.numberOfEntities(), rhs.entityLength()),
    data_(rhs.asEigenVector())
    {}

    void permute(const Eigen::PermutationMatrix<Eigen::Dynamic> &permutation) override {
        permuteMethod(permutation);
    }

    const EigenVecType& asEigenVector() const {
        return data_;
    }

    EigenVecType& asEigenVector() {
        return data_;
    }

    RefEigenVecType dataRef(long i) {
        return data_.segment(calculateIndex(i),entityLength());
    };

protected:
    bool operator==(const DataVector<Scalar> &other) const {
        return (data_ == other.data_);
    }

    bool operator!=(const DataVector<Scalar> &other) const {
        return !(*this == other);
    }

    void permuteMethod(const Eigen::PermutationMatrix<Eigen::Dynamic> &permutation) {
        assert(permutation.indices().size() == numberOfEntities()
               && "The permutation vector length must be equal to the number of entities");
        if(entityLength() > 1) {
            data_ =  adaptedToEntityLength(permutation) * data_;
        }
        else
            data_ = permutation* data_;
    }

    Eigen::PermutationMatrix<Eigen::Dynamic> adaptedToEntityLength(
            const Eigen::PermutationMatrix<Eigen::Dynamic> &permutation){

        Eigen::VectorXi raw(permutation.indices().size()*entityLength());

        for (int i = 0; i < permutation.indices().size(); ++i) {
            auto originIdx = i*entityLength();
            auto targetIdx = permutation.indices()[i]*entityLength();

            raw[originIdx+0] = targetIdx+0;
            raw[originIdx+1] = targetIdx+1;
            raw[originIdx+2] = targetIdx+2;
        }
        return Eigen::PermutationMatrix<Eigen::Dynamic>(raw);
    }

    EigenVecType data_;
};

#endif //INPSIGHTS_DATAVECTOR_H
