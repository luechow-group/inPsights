//
// Created by Michael Heuer on 29.10.17.
//

#include "SpinTypesVector.h"

using namespace Eigen;

SpinTypesVector::SpinTypesVector(long size)
        : AbstractVector(size),
          spinTypes_(VectorXi::Constant(size,int(Spin::SpinType::none)))
{}

SpinTypesVector::SpinTypesVector(const VectorXi& spinTypes)
        : AbstractVector(spinTypes.size()),
          spinTypes_(spinTypes)
{
    assert(spinTypes_.maxCoeff() <= int(Spin::SpinType::alpha));
    assert(spinTypes_.minCoeff() >= int(Spin::SpinType::beta));
}

SpinTypesVector::SpinTypesVector(unsigned long numberOfAlphaElectrons, unsigned long numberOfBetaElectrons)
        : SpinTypesVector(0)
{
    for (unsigned long i = 0; i < numberOfAlphaElectrons+numberOfBetaElectrons; ++i) {
        Spin::SpinType spinType;
        if (i < numberOfAlphaElectrons) spinType = Spin::SpinType::alpha;
        else  spinType = Spin::SpinType::beta;

        this->append(spinType);
    }
}

Spin::SpinType SpinTypesVector::operator[](long i) const {
    return  Spin::SpinType(spinTypes_[calculateIndex(i)]);
}

void SpinTypesVector::insert(Spin::SpinType spinType, long i) {
    assert(i >= 0 && "The index must be positive.");
    assert(i <= numberOfEntities() && "The index must be smaller than the number of entities.");

    VectorXi before = spinTypes_.head(i);
    VectorXi after = spinTypes_.tail(numberOfEntities()-i);

    spinTypes_.resize(numberOfEntities()+1);
    spinTypes_ << before, int(spinType), after;

    incrementNumberOfEntities();
}

void SpinTypesVector::prepend(Spin::SpinType spinType) {
    this->insert(spinType,0);
}

void SpinTypesVector::append(Spin::SpinType spinType) {
    this->insert(spinType,numberOfEntities());
}

const VectorXi& SpinTypesVector::spinTypesAsEigenVector() const {
    return spinTypes_;
}

VectorXi& SpinTypesVector::spinTypesAsEigenVector() {
    return spinTypes_;
}

void SpinTypesVector::permute(long i, long j) {
    if(i != j) {
        int temp = spinTypes_[calculateIndex(i)];
        spinTypes_[calculateIndex(i)] = spinTypes_[calculateIndex(j)];
        spinTypes_[calculateIndex(j)] = temp;
    }
}

std::ostream& operator<<(std::ostream& os, const SpinTypesVector& sc){
    for (unsigned long i = 0; i < sc.numberOfEntities(); i++) {
        os << Spin::toString(sc[i]) << std::endl;
    }
    return os;
}
