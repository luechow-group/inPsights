//
// Created by Michael Heuer on 30.10.17.
//

#include "ElectronCollections.h"

ElectronCollections::ElectronCollections()
        : ParticleCollections(),
          spinTypeCollection_(SpinTypeCollection())
{}

ElectronCollections::ElectronCollections(const std::vector<ElectronCollection> &electronCollections)
        : ParticleCollections(),
          spinTypeCollection_(SpinTypeCollection((*electronCollections.begin()).spinTypesAsEigenVector()))
{}

ElectronCollections::ElectronCollections(const std::vector<ParticleCollection> &particleCollections_,
                                         const SpinTypeCollection &spinTypeCollection)
        : ParticleCollections(particleCollections_),
          spinTypeCollection_(spinTypeCollection)
{}
