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
    explicit ElectronCollections(const SpinTypeCollection& spinTypeCollection);
    explicit ElectronCollections(const ElectronCollection& electronCollection);
    explicit ElectronCollections(const std::vector<ElectronCollection>& electronCollectionVector);
    explicit ElectronCollections(const PositionCollections& electronCollection);

    explicit ElectronCollections(const PositionCollections& electronCollection,
                                 const SpinTypeCollection& spinTypeCollection);

    ElectronCollection operator[](long i) const;

    const SpinTypeCollection& spinTypeCollection() const;
    SpinTypeCollection& spinTypeCollection();

    void insert (const ElectronCollection& electronCollection, long i);
    void append (const ElectronCollection& electronCollection);
    void prepend(const ElectronCollection& electronCollection);
    void permute(long i, long j) override;

private:
    SpinTypeCollection spinTypeCollection_;
};

#endif //AMOLQCGUI_ELECTRONCOLLECTIONPATH_H
