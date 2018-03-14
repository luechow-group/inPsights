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

double ParticlesVectorCollection::norm(long i, long j) const{
    return positionsVectorCollection_.norm(i, j);
}

