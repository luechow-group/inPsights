//
// Created by Michael Heuer on 29.10.17.
//

#include "SpinTypeCollection.h"

using namespace Eigen;

SpinTypeCollection::SpinTypeCollection(long size)
        :size_(size),
         spinTypes_(VectorXi::Constant(size,int(Spin::SpinType::none)))
{}

SpinTypeCollection::SpinTypeCollection(const VectorXi& spinTypes)
        : size_(spinTypes.size()),
          spinTypes_(size_)
{
    assert(spinTypes.minCoeff() >= int(Spin::SpinType::alpha));
    assert(spinTypes.maxCoeff() <= int(Spin::SpinType::none));

    spinTypes_ = spinTypes;
}

Spin::SpinType SpinTypeCollection::spinType(long i) {
    return  Spin::SpinType(spinTypes_[i]);
}

void SpinTypeCollection::setSpinType(long i, Spin::SpinType spinType) {
    spinTypes_[i] = int(spinType);
}

VectorXi SpinTypeCollection::asVectorXi() {
    return spinTypes_;
}
