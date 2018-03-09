//
// Created by Michael Heuer on 30.10.17.
//

#include "ParticlesVectorCollection.h"

ParticlesVectorCollection::ParticlesVectorCollection(const PositionsVectorCollection& positionsVectorCollection)
        : AbstractVector(positionsVectorCollection.numberOfEntities()),
          positionsVectorCollection_(positionsVectorCollection) {}

const PositionsVectorCollection &ParticlesVectorCollection::positionsVectorCollection() const {
    return positionsVectorCollection_;
}

PositionsVectorCollection &ParticlesVectorCollection::positionsVectorCollection() {
    return positionsVectorCollection_;
}
