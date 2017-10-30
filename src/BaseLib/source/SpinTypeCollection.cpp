//
// Created by Michael Heuer on 29.10.17.
//

#include <Particle.h>
#include "SpinTypeCollection.h"

using namespace Eigen;

SpinTypeCollection::SpinTypeCollection(long size)
        : numberOfSpinTypes_(size),
          spinTypes_(VectorXi::Constant(size,int(Spin::SpinType::none)))
{}

SpinTypeCollection::SpinTypeCollection(const VectorXi& spinTypes)
        : numberOfSpinTypes_(spinTypes.size()),
          spinTypes_(numberOfSpinTypes_)
{
    assert(spinTypes.minCoeff() >= int(Spin::SpinType::alpha));
    assert(spinTypes.maxCoeff() <= int(Spin::SpinType::none));

    spinTypes_ = spinTypes;
}

Spin::SpinType SpinTypeCollection::spinType(long i) {
    return  Spin::SpinType(spinTypes_[i]);
}

unsigned long SpinTypeCollection::numberOfSpinTypes() const {
    return numberOfSpinTypes_;
}

void SpinTypeCollection::insert(Spin::SpinType spinType, long i) {
    VectorXi before = spinTypes_.head(i);
    VectorXi after = spinTypes_.tail(numberOfSpinTypes_-i);

    spinTypes_.resize(numberOfSpinTypes_+1);
    spinTypes_ << before, int(spinType), after;
    ++numberOfSpinTypes_;
}

void SpinTypeCollection::prepend(Spin::SpinType spinType) {
    this->insert(spinType,0);
}

void SpinTypeCollection::append(Spin::SpinType spinType) {
    this->insert(spinType,numberOfSpinTypes_);
}

void SpinTypeCollection::setSpinType(long i, Spin::SpinType spinType) {
    spinTypes_[i] = int(spinType);
}

VectorXi SpinTypeCollection::spinTypesAsEigenVector() {
    return spinTypes_;
}
