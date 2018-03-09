
#include "AtomCollection.h"
#include "ToString.h"

using namespace Eigen;

AtomCollection::AtomCollection(const VectorXd &positions)
        : ParticleCollection(PositionsVector(positions)),
          elementTypeCollection_(numberOfEntities())
{}



AtomCollection::AtomCollection(const VectorXd &positions, const VectorXi &elementTypes)
        : ParticleCollection(PositionsVector(positions)),
          elementTypeCollection_(elementTypes) {

    assert(numberOfEntities() == positionsVector_.numberOfEntities()
           && numberOfEntities() == elementTypeCollection_.numberOfEntities()
           && "The number of entities in ParticleCollection, PositionsVector, and SpinTypeCollection must match.");
}

AtomCollection::AtomCollection(const PositionsVector &positionsVector,
                               const ElementTypeCollection &elementTypeCollection)
        : ParticleCollection(positionsVector),
          elementTypeCollection_(elementTypeCollection)
{}

Atom AtomCollection::operator[](long i) const {
    return Atom{positionsVector_[i], elementTypeCollection_[i]};
}

void AtomCollection::insert(const Atom& atom, long i) {
    assert(i >= 0 && "The index must be positive.");
    assert(i <= numberOfEntities() && "The index must be smaller than the number of entities.");

    positionsVector_.insert(atom.position(),i);
    elementTypeCollection_.insert(atom.elementType(),i);
    incrementNumberOfEntities();

    assert(positionsVector_.numberOfEntities() == numberOfEntities());
    assert(elementTypeCollection_.numberOfEntities() == numberOfEntities());
}

void AtomCollection::prepend(const Atom& atom) {
    this->insert(atom,0);
}

void AtomCollection::append(const Atom& atom) {
    this->insert(atom, numberOfEntities());
}

void AtomCollection::permute(long i, long j) {
    positionsVector_.permute(i,j);
    elementTypeCollection_.permute(i,j);
}

const ElementTypeCollection &AtomCollection::elementTypeCollection() const {
    return elementTypeCollection_;
}

ElementTypeCollection &AtomCollection::elementTypeCollection() {
    return elementTypeCollection_;
}

std::ostream& operator<<(std::ostream& os, const AtomCollection& ac){
    for (unsigned long i = 0; i < ac.numberOfEntities(); i++) {
        os << ToString::intToString(i + 1) << " " << ac[i] << std::endl;
    }
    return os;
}