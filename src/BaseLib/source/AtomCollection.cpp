
#include "AtomCollection.h"
#include "ToString.h"

using namespace Eigen;

AtomCollection::AtomCollection(const VectorXd &positions)
        : ParticlesVector(PositionsVector(positions)),
          elementTypesVector_(numberOfEntities())
{}



AtomCollection::AtomCollection(const VectorXd &positions, const VectorXi &elementTypes)
        : ParticlesVector(PositionsVector(positions)),
          elementTypesVector_(elementTypes) {

    assert(numberOfEntities() == positionsVector_.numberOfEntities()
           && numberOfEntities() == elementTypesVector_.numberOfEntities()
           && "The number of entities in ParticlesVector, PositionsVector, and SpinTypesVector must match.");
}

AtomCollection::AtomCollection(const PositionsVector &positionsVector,
                               const ElementTypesVector &elementTypesVector)
        : ParticlesVector(positionsVector),
          elementTypesVector_(elementTypesVector)
{}

Atom AtomCollection::operator[](long i) const {
    return Atom{positionsVector_[i], elementTypesVector_[i]};
}

void AtomCollection::insert(const Atom& atom, long i) {
    assert(i >= 0 && "The index must be positive.");
    assert(i <= numberOfEntities() && "The index must be smaller than the number of entities.");

    positionsVector_.insert(atom.position(),i);
    elementTypesVector_.insert(atom.elementType(),i);
    incrementNumberOfEntities();

    assert(positionsVector_.numberOfEntities() == numberOfEntities());
    assert(elementTypesVector_.numberOfEntities() == numberOfEntities());
}

void AtomCollection::prepend(const Atom& atom) {
    this->insert(atom,0);
}

void AtomCollection::append(const Atom& atom) {
    this->insert(atom, numberOfEntities());
}

void AtomCollection::permute(long i, long j) {
    positionsVector_.permute(i,j);
    elementTypesVector_.permute(i,j);
}

const ElementTypesVector &AtomCollection::elementTypesVector() const {
    return elementTypesVector_;
}

ElementTypesVector &AtomCollection::elementTypesVector() {
    return elementTypesVector_;
}

std::ostream& operator<<(std::ostream& os, const AtomCollection& ac){
    for (unsigned long i = 0; i < ac.numberOfEntities(); i++) {
        os << ToString::intToString(i + 1) << " " << ac[i] << std::endl;
    }
    return os;
}