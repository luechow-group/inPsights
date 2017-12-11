//
// Created by Michael Heuer on 29.10.17.
//

#include "ElectronCollection.h"

ElectronCollection::ElectronCollection(const Eigen::VectorXd &positions)
        : ParticleCollection(positions),
          SpinTypeCollection(this->numberOfParticles()) {}

ElectronCollection::ElectronCollection(const Eigen::VectorXd &positions, const Eigen::VectorXi &spinTypes)
        : ParticleCollection(positions),
          SpinTypeCollection(spinTypes) {

    assert(this->numberOfParticles() == spinTypes.size()
           && "The number of particles in ParticleCollection and the number of spin type in SpinTypeCollection must match.");
}

/*
ElectronCollection::ElectronCollection(const ParticleCollection &particleCollection,
                                       const SpinTypeCollection &spinTypeCollection)
        : ElectronCollection(particleCollection.positionsAsEigenVector(),
                             spinTypeCollection.spinTypesAsEigenVector()){
}*/

ElectronCollection::ElectronCollection(const ParticleCollection &particleCollection,
                                       const SpinTypeCollection &spinTypeCollection)
        : ParticleCollection(particleCollection),
          SpinTypeCollection(spinTypeCollection)
{}

Electron ElectronCollection::electron(long i) {
    Particle particle = (*this)[i];
    return Electron(particle, spinType(i));
}

void ElectronCollection::insert(const Electron& electron, long i) {
    ParticleCollection::insert(static_cast<Particle>(electron),i);
    SpinTypeCollection::insert(electron.spinType(),i);
    assert(numberOfParticles() == numberOfSpinTypes());
}

void ElectronCollection::prepend(const Electron& electron) {
    this->insert(electron,0);
}

void ElectronCollection::append(const Electron& electron) {
    this->insert(electron, ParticleCollection::numberOfParticles_);
}

void ElectronCollection::permute(long i, long j) {
    ParticleCollection::permute(i,j);
    SpinTypeCollection::permute(i,j);
}
