//
// Created by Michael Heuer on 29.10.17.
//

#ifndef AMOLQCGUI_PARTICLECOLLECTION_H
#define AMOLQCGUI_PARTICLECOLLECTION_H

#include "AbstractVector.h"
#include "Particle.h"
#include "PositionsVector.h"

/*
 * ParticleCollection only exists as an abstract interface to more specialized collections.
 */
class ParticleCollection : public AbstractVector{
public:
    Particle particle(long i) const;

    const PositionsVector & positionsVector() const;
    PositionsVector & positionsVector();

protected:
    PositionsVector positionsVector_;

    void permute(long i, long j) override = 0;
    ParticleCollection() = default;
    explicit ParticleCollection(const PositionsVector& positionsVector);
};

#endif //AMOLQCGUI_PARTICLECOLLECTION_H
