//
// Created by Michael Heuer on 30.10.17.
//

#include "ParticleCollections.h"

ParticleCollections::ParticleCollections()
        : numberOfParticleColledctions_(0),
          particleCollections_()
{}

ParticleCollections::ParticleCollections(std::vector<ParticleCollection> particleCollections)
        : numberOfParticleColledctions_(particleCollections.size()),
          particleCollections_(particleCollections)
{}

ParticleCollection ParticleCollections::operator[](long i) {
    return particleCollections_[i];
}

void ParticleCollections::insert(const ParticleCollection &particleCollection, long i) {
    particleCollections_.insert(particleCollections_.begin()+i,particleCollection);
}

void ParticleCollections::append(const ParticleCollection &particleCollection) {
    //this->insert(particleCollection,0);
    particleCollections_.push_back(particleCollection);
}

void ParticleCollections::prepend(const ParticleCollection &particleCollection) {
    this->insert(particleCollection, particleCollection.numberOfParticles() );
}

std::vector<ParticleCollection> ParticleCollections::getParticleCollections() const {
    return particleCollections_;
}
