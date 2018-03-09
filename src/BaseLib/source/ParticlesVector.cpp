//
// Created by Michael Heuer on 29.10.17.
//

#include "ParticlesVector.h"

ParticlesVector::ParticlesVector(const PositionsVector &positionsVector)
        : AbstractVector(positionsVector.numberOfEntities()),
          positionsVector_(positionsVector){}

Particle ParticlesVector::particle(long i) const {
    return Particle(positionsVector_[i]);
}

const PositionsVector & ParticlesVector::positionsVector() const {
    return positionsVector_;
}

PositionsVector & ParticlesVector::positionsVector(){
    return positionsVector_;
}
