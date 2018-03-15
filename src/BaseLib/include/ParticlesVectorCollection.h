//
// Created by Michael Heuer on 30.10.17.
//

#ifndef AMOLQCPP_PARTICLECOLLECTIONS_H
#define AMOLQCPP_PARTICLECOLLECTIONS_H

#include <vector>
#include "ParticlesVector.h"
#include "PositionsVectorCollection.h"

class ParticlesVectorCollection : public AbstractVector{
public:
    const PositionsVectorCollection& positionsVectorCollection() const;
    PositionsVectorCollection& positionsVectorCollection();
    double norm(long i, long j) const;

protected:
    PositionsVectorCollection positionsVectorCollection_;

    void permute(long i, long j) override = 0;
    ParticlesVectorCollection() = default;
    explicit ParticlesVectorCollection(const PositionsVectorCollection& positionsVectorCollection);

//private:
//    long calculateIndex(long i) const final = 0;
};

#endif //AMOLQCPP_PARTICLECOLLECTIONS_H
