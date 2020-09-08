// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_INSERTABLEVECTOR_H
#define INPSIGHTS_INSERTABLEVECTOR_H

#include "DataVector.h"

template<typename Scalar>
class InsertableVector : public DataVector<Scalar> {
protected:
    using EigenVecType = Eigen::Matrix<Scalar, Eigen::Dynamic,1>;

    InsertableVector(long numberOfEntities = 0, long entityLength = 1)
    :
    DataVector<Scalar>(numberOfEntities,entityLength) {}

    void insert(const EigenVecType& data, long i) { // For types use vector VectorXi
        assert(i >= 0 && "The index must be positive.");
        assert(i <= AbstractVector::numberOfEntities() && "The index must be smaller than the number of entities.");

        long start = AbstractVector::calculateIndex(i);

        EigenVecType before = DataVector<Scalar>::data_.head(start);
        EigenVecType after = DataVector<Scalar>::data_.tail(AbstractVector::numberOfEntities()*AbstractVector::entityLength()-start);

        DataVector<Scalar>::data_.resize(AbstractVector::numberOfEntities()*AbstractVector::entityLength()+AbstractVector::entityLength());
        DataVector<Scalar>::data_ << before, data, after;

        AbstractVector::incrementNumberOfEntities();
    }
};

#endif //INPSIGHTS_INSERTABLEVECTOR_H
