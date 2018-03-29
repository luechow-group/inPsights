//
// Created by Michael Heuer on 29.10.17.
//

#ifndef AMOLQCPP_PARTICLESVECTOR_H
#define AMOLQCPP_PARTICLESVECTOR_H

#include "AbstractVector.h"
#include "Particle.h"
#include "PositionsVector.h"

/*
 * ParticlesVector only exists as an abstract interface to more specialized collections.
 */
class ParticlesVector : public AbstractVector{
public:
    Particle particle(long i) const;

    const PositionsVector & positionsVector() const;
    PositionsVector & positionsVector();

protected:
    PositionsVector positionsVector_;

    void permute(long i, long j) override = 0;
    ParticlesVector() = default;
    explicit ParticlesVector(const PositionsVector& positionsVector);
};

#endif //AMOLQCPP_PARTICLESVECTOR_H
