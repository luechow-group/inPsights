
#include "AtomCollection.h"

using namespace Eigen;

AtomCollection::AtomCollection(const VectorXd &positions)
        : ParticleCollection(positions),
          ElementTypeCollection(ParticleCollection::numberOfParticles())
{}



AtomCollection::AtomCollection(const VectorXd &positions, const VectorXi &elementTypes)
        : ParticleCollection(positions),
          ElementTypeCollection(elementTypes)
{}

AtomCollection::AtomCollection(const ParticleCollection& particleCollection,
                               const ElementTypeCollection& elementTypeCollection)
        : ParticleCollection(particleCollection),
          ElementTypeCollection(elementTypeCollection)
{}

Atom AtomCollection::atom(long i) {
    Particle particle = (*this)[i];
    return Atom(particle, elementType(i));
}

void AtomCollection::insert(const Atom& atom, long i) {
    ParticleCollection::insert(static_cast<Particle>(atom),i);
    ElementTypeCollection::insert(atom.elementType(),i);
    assert(numberOfParticles() == numberOfElementTypes());
}

void AtomCollection::prepend(const Atom& atom) {
    this->insert(atom,0);
}

void AtomCollection::append(const Atom& atom) {
    this->insert(atom, ParticleCollection::numberOfParticles());
}


void AtomCollection::addAtom(double x, double y, double z, const Elements::ElementType &elementType) {
  append(Atom(Eigen::Vector3d(x,y,z), elementType));
};

void AtomCollection::addAtom(const Vector3d &position, const Elements::ElementType &elementType) {
  append(Atom(position, elementType));
};

