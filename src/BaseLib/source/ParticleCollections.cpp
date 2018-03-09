//
// Created by Michael Heuer on 30.10.17.
//

#include "ParticleCollections.h"

ParticleCollections::ParticleCollections(const PositionCollections& positionCollections)
        : AbstractCollection(positionCollections.numberOfEntities()),
          positionCollections_(positionCollections) {}

const PositionCollections &ParticleCollections::positionCollections() const {
    return positionCollections_;
}

PositionCollections &ParticleCollections::positionCollections() {
    return positionCollections_;
}
