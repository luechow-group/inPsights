//
// Created by Michael Heuer on 30.10.17.
//

#ifndef AMOLQCGUI_PARTICLECOLLECTIONS_H
#define AMOLQCGUI_PARTICLECOLLECTIONS_H

#include <vector>
#include "ParticleCollection.h"

class ParticleCollections{
public:
    ParticleCollections();
    ParticleCollections(std::vector<ParticleCollection> particleCollections);

    ParticleCollection operator[](long i);

    void insert (const ParticleCollection& particleCollection, long i);
    void append (const ParticleCollection& particleCollection);
    void prepend(const ParticleCollection& particleCollection);


    std::vector<ParticleCollection> getParticleCollections() const;

private:
    unsigned long numberOfParticleColledctions_;
    std::vector<ParticleCollection> particleCollections_;
};

#endif //AMOLQCGUI_PARTICLECOLLECTIONS_H
