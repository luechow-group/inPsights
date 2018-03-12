
#include "AtomsVector.h"
#include "ToString.h"

using namespace Eigen;

AtomsVector::AtomsVector(const VectorXd &positions)
        : ParticlesVector(PositionsVector(positions)),
          elementTypesVector_(numberOfEntities())
{}

AtomsVector::AtomsVector(const VectorXd &positions, const VectorXi &elementTypes)
        : ParticlesVector(PositionsVector(positions)),
          elementTypesVector_(elementTypes) {

    assert(numberOfEntities() == positionsVector_.numberOfEntities()
           && numberOfEntities() == elementTypesVector_.numberOfEntities()
           && "The number of entities in ParticlesVector, PositionsVector, and SpinTypesVector must match.");
}

AtomsVector::AtomsVector(const PositionsVector &positionsVector,
                               const ElementTypesVector &elementTypesVector)
        : ParticlesVector(positionsVector),
          elementTypesVector_(elementTypesVector)
{}

Atom AtomsVector::operator[](long i) const {
    return Atom{positionsVector_[i], elementTypesVector_[i]};
}

void AtomsVector::insert(const Atom& atom, long i) {
    assert(i >= 0 && "The index must be positive.");
    assert(i <= numberOfEntities() && "The index must be smaller than the number of entities.");

    positionsVector_.insert(atom.position(),i);
    elementTypesVector_.insert(atom.elementType(),i);
    incrementNumberOfEntities();

    assert(positionsVector_.numberOfEntities() == numberOfEntities());
    assert(elementTypesVector_.numberOfEntities() == numberOfEntities());
}

void AtomsVector::prepend(const Atom& atom) {
    this->insert(atom,0);
}

void AtomsVector::append(const Atom& atom) {
    this->insert(atom, numberOfEntities());
}

void AtomsVector::permute(long i, long j) {
    positionsVector_.permute(i,j);
    elementTypesVector_.permute(i,j);
}

const ElementTypesVector &AtomsVector::elementTypesVector() const {
    return elementTypesVector_;
}

ElementTypesVector &AtomsVector::elementTypesVector() {
    return elementTypesVector_;
}

std::ostream& operator<<(std::ostream& os, const AtomsVector& ac){
    for (unsigned long i = 0; i < ac.numberOfEntities(); i++) {
        os << ToString::unsignedLongToString(i + 1) << " " << ac[i] << std::endl;
    }
    return os;
}