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
    ParticleCollections(const std::vector<ParticleCollection> &particleCollections);

    ParticleCollection operator[](long i) const;

    void insert (const ParticleCollection& particleCollection, long i);
    void append (const ParticleCollection& particleCollection);
    void prepend(const ParticleCollection& particleCollection);

    //TODO remove: std::vector<ParticleCollection>::const_iterator begin() const { return particleCollections_.begin(); };

    std::vector<ParticleCollection> getParticleCollections() const;

    unsigned long getNumberOfParticleCollections() const;
    unsigned long getNumberOfParticles() const;

private:
    unsigned long numberOfParticles_;
    std::vector<ParticleCollection> particleCollections_;
};

#endif //AMOLQCGUI_PARTICLECOLLECTIONS_H
