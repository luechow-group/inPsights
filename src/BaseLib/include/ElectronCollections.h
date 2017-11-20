//
// Created by Michael Heuer on 30.10.17.
//

#ifndef AMOLQCGUI_ELECTRONCOLLECTIONPATH_H
#define AMOLQCGUI_ELECTRONCOLLECTIONPATH_H

#include "ParticleCollections.h"
#include "ElectronCollection.h"
#include "SpinTypeCollection.h"

class ElectronCollections : public ParticleCollections{
public:
    ElectronCollections();
    ElectronCollections(const SpinTypeCollection& spinTypeCollection);
    explicit ElectronCollections(const std::vector<ElectronCollection>& electronCollections);
    explicit ElectronCollections(const std::vector<ParticleCollection>& particleCollections);
    ElectronCollections(const std::vector<ParticleCollection>& particleCollections,
                        const SpinTypeCollection& spinTypeCollection);

    ElectronCollection getElectronCollection(long i) const;
    SpinTypeCollection getSpinTypeCollection() const;

    void insert (const ElectronCollection& electronCollection, long i);
    virtual void append (const ElectronCollection& electronCollection);
    void prepend(const ElectronCollection& electronCollection);

private:
    SpinTypeCollection spinTypeCollection_;
};

#endif //AMOLQCGUI_ELECTRONCOLLECTIONPATH_H
