//
// Created by Michael Heuer on 30.10.17.
//

#include <ElectronCollection.h>
#include "ParticleCollections.h"

ParticleCollections::ParticleCollections()
        : numberOfParticles_(0),
          particleCollections_()
{}

ParticleCollections::ParticleCollections(const std::vector<ParticleCollection> &particleCollections)
        : numberOfParticles_(0),
          particleCollections_(particleCollections)
{
    if( !particleCollections.empty()) {
        numberOfParticles_ = particleCollections[0].numberOfParticles();

        for (const auto &particleCollection : particleCollections) {
            assert( numberOfParticles_ == particleCollection.numberOfParticles()
                    && "All particle collections must contain the same number of particles.");
        }
    }
}

ParticleCollection ParticleCollections::operator[](long i) const {
    return particleCollections_[i];
}

ParticleCollection ParticleCollections::front() const {
    return particleCollections_.front();
}

ParticleCollection ParticleCollections::back() const {
    return particleCollections_.back();
}

void ParticleCollections::insert(const ParticleCollection &particleCollection, long i) {

    if(length() == 0) {
        assert( i == 0 && "If the collection is empty, the index i must be zero ");
        particleCollections_ = {particleCollection};
        numberOfParticles_ = particleCollection.numberOfEntities();
    } else {
        assert(numberOfParticles_ == particleCollection.numberOfParticles()
               && "All particle collections must contain the same number of particles.");
    }
    particleCollections_.insert(particleCollections_.begin()+i,particleCollection);

}

void ParticleCollections::append(const ParticleCollection &particleCollection) {
    if(length() == 0) {
        numberOfParticles_ = particleCollection.numberOfParticles();
    }
    else {
        assert( numberOfParticles_ == particleCollection.numberOfParticles()
                && "All particle collections must contain the same number of particles.");
    }
    particleCollections_.emplace_back(particleCollection);
}

void ParticleCollections::prepend(const ParticleCollection &particleCollection) {
    this->insert(particleCollection,0);
}

std::vector<ParticleCollection> ParticleCollections::getParticleCollections() const {
    return particleCollections_;
}

unsigned long ParticleCollections::length() const {
    return particleCollections_.size();
}

unsigned long ParticleCollections::getNumberOfParticles() const {
    return numberOfParticles_;
}
