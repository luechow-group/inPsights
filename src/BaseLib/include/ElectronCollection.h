//
// Created by Michael Heuer on 29.10.17.
//

#ifndef AMOLQCGUI_ELECTRONCOLLECTION_H
#define AMOLQCGUI_ELECTRONCOLLECTION_H

#include "ParticleCollection.h"
#include "SpinTypeCollection.h"
#include "Electron.h"

class ElectronCollection : public ParticleCollection, public SpinTypeCollection{
public:
    ElectronCollection() = default;
    explicit ElectronCollection(const VectorXd& positions);
    ElectronCollection(const VectorXd& positions, const VectorXi& spinTypes);

    ElectronCollection(const ParticleCollection& particleCollection,
                       const SpinTypeCollection& spinTypeCollection);

    Electron electron(long i);

    void insert (const Electron& electron, long i);
    void append (const Electron& electron);
    void prepend(const Electron& electron);
};

#endif //AMOLQCGUI_ELECTRONCOLLECTION_H
