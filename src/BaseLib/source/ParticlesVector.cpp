//
// Created by Michael Heuer on 29.10.17.
//p

#include "ParticlesVector.h"
#include "ToString.h"

ParticlesVector::ParticlesVector(const PositionsVector &positionsVector)
        : AbstractVector(positionsVector.numberOfEntities()),
          positionsVector_(positionsVector),
          typesVector_(numberOfEntities())
{}

ParticlesVector::ParticlesVector(const PositionsVector &positionsVector,
                                 const TypesVector &typesVector)
        : AbstractVector(positionsVector.numberOfEntities()),
          positionsVector_(positionsVector),
          typesVector_(typesVector) {
    assert(numberOfEntities() == positionsVector_.numberOfEntities()
           && numberOfEntities() == typesVector_.numberOfEntities()
           && "The number of entities in ParticlesVector, PositionsVector, and TypesVector must match.");
}

Particle ParticlesVector::particle(long i) const {
    return {positionsVector_[i],typesVector_.type(i)};
}

const PositionsVector & ParticlesVector::positionsVector() const {
    return positionsVector_;
}

PositionsVector & ParticlesVector::positionsVector(){
    return positionsVector_;
}

const TypesVector & ParticlesVector::typesVector() const {
    return typesVector_;
}

TypesVector & ParticlesVector::typesVector(){
    return typesVector_;
}

void ParticlesVector::insert(const Particle& particle, long i) {
    assert(i >= 0 && "The index must be positive.");
    assert(i <= numberOfEntities() && "The index must be smaller than the number of entities.");

    positionsVector_.insert(particle.position(),i);
    typesVector_.insert(particle.type(),i);
    incrementNumberOfEntities();

    assert(positionsVector_.numberOfEntities() == numberOfEntities());
    assert(typesVector_.numberOfEntities() == numberOfEntities());
}

void ParticlesVector::prepend(const Particle& particle) {
    this->insert(particle,0);
}

void ParticlesVector::append(const Particle& particle) {
    this->insert(particle,numberOfEntities());
}

void ParticlesVector::permute(long i, long j) {
    positionsVector_.permute(i,j);
    typesVector_.permute(i,j);
}

std::ostream& operator<<(std::ostream& os, const ParticlesVector& pv){
    for (unsigned long i = 0; i < pv.numberOfEntities(); i++) {
        os << ToString::unsignedLongToString(i + 1) << " " << pv.particle(i) << std::endl;
    }
    return os;
}
