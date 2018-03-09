//
// Created by Michael Heuer on 30.10.17.
//

#ifndef AMOLQCGUI_PARTICLECOLLECTIONS_H
#define AMOLQCGUI_PARTICLECOLLECTIONS_H

#include <vector>
#include "ParticleCollection.h"
#include "PositionCollections.h"

class ParticleCollections : public AbstractCollection{
public:
    const PositionCollections& positionCollections() const;
    PositionCollections& positionCollections();

protected:
    PositionCollections positionCollections_;

    void permute(long i, long j) override = 0;
    ParticleCollections() = default;
    explicit ParticleCollections(const PositionCollections& positionCollections);

//private:
//    long calculateIndex(long i) const final = 0;
};

#endif //AMOLQCGUI_PARTICLECOLLECTIONS_H
