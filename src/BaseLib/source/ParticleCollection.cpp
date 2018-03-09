//
// Created by Michael Heuer on 29.10.17.
//

#include "ParticleCollection.h"

ParticleCollection::ParticleCollection(const PositionsVector &positionsVector)
        : AbstractCollection(positionsVector.numberOfEntities()),
          positionsVector_(positionsVector){}

Particle ParticleCollection::particle(long i) const {
    return Particle(positionsVector_[i]);
}

const PositionsVector & ParticleCollection::positionsVector() const {
    return positionsVector_;
}

PositionsVector & ParticleCollection::positionsVector(){
    return positionsVector_;
}
