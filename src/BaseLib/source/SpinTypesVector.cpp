//
// Created by Michael Heuer on 29.10.17.
//

#include "ParticlesVector.h"

using namespace Eigen;

SpinTypesVector::SpinTypesVector(long size)
        : TypesVector(size)
{}

SpinTypesVector::SpinTypesVector(const VectorXi& types)
        : TypesVector(types)
{
    assert(types_.maxCoeff() <= Spin::spinTypeToInt(Spin::SpinType::alpha));
    assert(types_.minCoeff() >= Spin::spinTypeToInt(Spin::SpinType::beta));
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
    return Spin::spinTypeFromInt(TypesVector::type(i));
}

void SpinTypesVector::insert(Spin::SpinType spinType, long i) {
    TypesVector::insert(Spin::spinTypeToInt(spinType),i);
}

void SpinTypesVector::prepend(Spin::SpinType spinType) {
    TypesVector::prepend(Spin::spinTypeToInt(spinType));
}

void SpinTypesVector::append(Spin::SpinType spinType) {
    TypesVector::append(Spin::spinTypeToInt(spinType));
}

std::ostream& operator<<(std::ostream& os, const SpinTypesVector& sc){
    for (unsigned long i = 0; i < sc.numberOfEntities(); i++) {
        os << Spin::toString(sc[i]) << std::endl;
    }
    return os;
}
