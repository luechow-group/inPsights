//
// Created by Michael Heuer on 29.10.17.
//

#ifndef AMOLQCGUI_PARTICLECOLLECTION_H
#define AMOLQCGUI_PARTICLECOLLECTION_H

#include "AbstractCollection.h"
#include "Particle.h"
#include "PositionCollection.h"


/*
 * ParticleCollection only exists as an abstract interface to more specialized collections.
 */
class ParticleCollection : public AbstractCollection{
public:
    Particle particle(long i) const;

    const PositionCollection & positionCollection() const;
    PositionCollection & positionCollection();

protected:
    PositionCollection positionCollection_;

    void permute(long i, long j) override = 0;
    ParticleCollection() = default;
    explicit ParticleCollection(const PositionCollection& positionCollection);
};

#endif //AMOLQCGUI_PARTICLECOLLECTION_H
