
#include "AtomCollection.h"

using namespace Eigen;

AtomCollection::AtomCollection(const VectorXd &positions)
        : ParticleCollection(positions),
          ElementTypeCollection(ParticleCollection::numberOfParticles())
{}

AtomCollection::AtomCollection(const VectorXd &positions, const VectorXi &spinTypes)
        : ParticleCollection(positions),
          ElementTypeCollection(spinTypes)
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
    this->insert(atom, ParticleCollection::numberOfParticles_);
}


void AtomCollection::addAtom(const Elements::ElementType &elementType,
                             const double x, const double y, const double z) {
  append(Atom(Eigen::Vector3d(x,y,z), elementType));
};

void AtomCollection::addAtom(const Elements::ElementType &elementType,
                             const Eigen::Vector3d &position) {
  append(Atom(position, elementType));
};


