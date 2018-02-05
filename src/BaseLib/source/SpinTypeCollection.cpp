//
// Created by Michael Heuer on 29.10.17.
//

#include <Particle.h>
#include "SpinTypeCollection.h"
#include "SpinType.h"

using namespace Eigen;

SpinTypeCollection::SpinTypeCollection(unsigned long size)
        : numberOfSpinTypes_(size),
          spinTypes_(VectorXi::Constant(size,int(Spin::SpinType::none)))
{}

SpinTypeCollection::SpinTypeCollection(const VectorXi& spinTypes)
        : numberOfSpinTypes_(spinTypes.size()),
          spinTypes_(numberOfSpinTypes_)
{
    assert(spinTypes.maxCoeff() <= int(Spin::SpinType::alpha));
    assert(spinTypes.minCoeff() >= int(Spin::SpinType::beta));

    spinTypes_ = spinTypes;
}

SpinTypeCollection::SpinTypeCollection(unsigned long numberOfAlphaElectrons, unsigned long numberOfBetaElectrons)
        : SpinTypeCollection(0)
{
    for (unsigned long i = 0; i < numberOfAlphaElectrons+numberOfBetaElectrons; ++i) {
        Spin::SpinType spinType;
        if (i < numberOfAlphaElectrons) spinType = Spin::SpinType::alpha;
        else  spinType = Spin::SpinType::beta;

        this->append(spinType);
    }
}

Spin::SpinType SpinTypeCollection::spinType(long i) const {
    return  Spin::SpinType(spinTypes_[i]);
}

unsigned long SpinTypeCollection::numberOfSpinTypes() const {
    return numberOfSpinTypes_;
}

void SpinTypeCollection::insert(Spin::SpinType spinType, long i) {
    VectorXi before = spinTypes_.head(i);
    VectorXi after = spinTypes_.tail(numberOfSpinTypes_-i);

    spinTypes_.resize(numberOfSpinTypes_+1);
    //spinTypes_ << before, int(spinType), after;
    spinTypes_.head(i) = before;
    spinTypes_.segment(i,1) = Eigen::Matrix<int,1,1>(int(spinType));
    spinTypes_.tail(numberOfSpinTypes_-i) = after;
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

VectorXi SpinTypeCollection::spinTypesAsEigenVector() const {
    return spinTypes_;
}

void SpinTypeCollection::permute(long i, long j) {
    assert( i >= 0 && i < numberOfSpinTypes_
            && "Index i must be greater than zero and smaller than the number of spins." );
    assert( j >= 0 && j < numberOfSpinTypes_
            && "Index j must be greater than zero and smaller than the number of spins." );
    if(i != j) {
        int temp = spinTypes_[i];
        spinTypes_[i] = spinTypes_[j];
        spinTypes_[j] = temp;
    }
}

std::ostream& operator<<(std::ostream& os, const SpinTypeCollection& sc){
    for (unsigned long i = 0; i < sc.numberOfSpinTypes(); i++) {
        os << std::string(ParticleFormat::significantDigits+1, ' ')  << Spin::toString(sc.spinType(i)) << ParticleFormat::separator;
    }
    std::cout << std::endl;
    return os;
}
