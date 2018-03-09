//
// Created by Michael Heuer on 29.10.17.
//

#include "ParticleCollection.h"

ParticleCollection::ParticleCollection(const PositionCollection &positionCollection)
        : AbstractCollection(positionCollection.numberOfEntities()),
          positionCollection_(positionCollection){}

Particle ParticleCollection::particle(long i) const {
    return Particle(positionCollection_[i]);
}

const PositionCollection & ParticleCollection::positionCollection() const {
    return positionCollection_;
}

PositionCollection & ParticleCollection::positionCollection(){
    return positionCollection_;
}
