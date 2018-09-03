//
// Created by Michael Heuer on 03.09.18.
//

#ifndef AMOLQCPP_IINSERTABLE_H
#define AMOLQCPP_IINSERTABLE_H

#include "ISliceable.h"

template<typename Scalar>
class IInsertable : public ISliceable<Scalar> {
public:

    IInsertable(long numberOfEntities = 0, long entityLength = 1)
            :
            ISliceable<Scalar>(numberOfEntities,entityLength) {}

    using VectorType = Eigen::Matrix<Scalar, Eigen::Dynamic,1>;

protected:
    void insert(const VectorType& data, long i) { // For types use vector VectorXi
        assert(i >= 0 && "The index must be positive.");
        assert(i <= AbstractVector::numberOfEntities() && "The index must be smaller than the number of entities.");

        long start = AbstractVector::calculateIndex(i);

        VectorType before = ISliceable<Scalar>::data_.head(start);
        VectorType after = ISliceable<Scalar>::data_.tail(AbstractVector::numberOfEntities()*AbstractVector::entityLength()-start);

        ISliceable<Scalar>::data_.resize(AbstractVector::numberOfEntities()*AbstractVector::entityLength()+AbstractVector::entityLength());
        ISliceable<Scalar>::data_ << before, data, after;

        AbstractVector::incrementNumberOfEntities();
        ISliceable<Scalar>::resetRef(); // reset, because slices and refs are invalid now
    }
};

#endif //AMOLQCPP_IINSERTABLE_H
