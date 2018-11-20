//
// Created by Michael Heuer on 03.09.18.
//

#ifndef INPSIGHTS_INSERTABLEVECTOR_H
#define INPSIGHTS_INSERTABLEVECTOR_H

#include "SliceableDataVector.h"

template<typename Scalar>
class InsertableVector : public SliceableDataVector<Scalar> {
protected:
    using EigenVecType = Eigen::Matrix<Scalar, Eigen::Dynamic,1>;

    InsertableVector(long numberOfEntities = 0, long entityLength = 1)
    :
    SliceableDataVector<Scalar>(numberOfEntities,entityLength) {}

    void insert(const EigenVecType& data, long i) { // For types use vector VectorXi
        assert(i >= 0 && "The index must be positive.");
        assert(i <= AbstractVector::numberOfEntities() && "The index must be smaller than the number of entities.");

        long start = AbstractVector::calculateIndex(i);

        EigenVecType before = SliceableDataVector<Scalar>::data_.head(start);
        EigenVecType after = SliceableDataVector<Scalar>::data_.tail(AbstractVector::numberOfEntities()*AbstractVector::entityLength()-start);

        SliceableDataVector<Scalar>::data_.resize(AbstractVector::numberOfEntities()*AbstractVector::entityLength()+AbstractVector::entityLength());
        SliceableDataVector<Scalar>::data_ << before, data, after;

        AbstractVector::incrementNumberOfEntities();
        SliceableDataVector<Scalar>::resetRef(); // reset, because slices and refs are invalid now
    }
};

#endif //INPSIGHTS_INSERTABLEVECTOR_H
