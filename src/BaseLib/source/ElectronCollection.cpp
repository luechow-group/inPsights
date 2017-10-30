//
// Created by Michael Heuer on 29.10.17.
//

#include "ElectronCollection.h"

ElectronCollection::ElectronCollection(const VectorXd &positions)
        : ParticleCollection(positions),
          SpinTypeCollection(this->numberOfParticles())
{}

ElectronCollection::ElectronCollection(const VectorXd &positions, const VectorXi &spinTypes)
        : ParticleCollection(positions),
          SpinTypeCollection(spinTypes)
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