//
// Created by Michael Heuer on 29.10.17.
//

#include "ElectronCollection.h"

ElectronCollection::ElectronCollection(const VectorXd &positions)
        : ParticleCollection(positions),
          SpinTypeCollection(ParticleCollection::size())
{}

ElectronCollection::ElectronCollection(const VectorXd &positions, const VectorXi &spinTypes)
        : ParticleCollection(positions),
          SpinTypeCollection(spinTypes)
{}

Electron ElectronCollection::electron(long i) {
    Particle particle = (*this)[i];
    return Electron(particle, spinType(i));
}
