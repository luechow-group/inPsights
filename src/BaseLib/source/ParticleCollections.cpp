//
// Created by Michael Heuer on 30.10.17.
//

#include "ParticleCollections.h"

ParticleCollections::ParticleCollections(const PositionsVectorCollection& positionsVectorCollection)
        : AbstractCollection(positionsVectorCollection.numberOfEntities()),
          positionsVectorCollection_(positionsVectorCollection) {}

const PositionsVectorCollection &ParticleCollections::positionsVectorCollection() const {
    return positionsVectorCollection_;
}

PositionsVectorCollection &ParticleCollections::positionsVectorCollection() {
    return positionsVectorCollection_;
}
