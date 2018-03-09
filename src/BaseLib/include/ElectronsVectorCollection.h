//
// Created by Michael Heuer on 30.10.17.
//

#ifndef AMOLQCGUI_ELECTRONCOLLECTIONPATH_H
#define AMOLQCGUI_ELECTRONCOLLECTIONPATH_H

#include "ParticleCollections.h"
#include "ElectronCollection.h"
#include "SpinTypeCollection.h"

class ElectronsVectorCollection : public ParticleCollections{
public:
    ElectronsVectorCollection();
    explicit ElectronsVectorCollection(const SpinTypeCollection& spinTypeCollection);
    explicit ElectronsVectorCollection(const ElectronCollection& electronCollection);
    explicit ElectronsVectorCollection(const std::vector<ElectronCollection>& electronCollectionVector);
    explicit ElectronsVectorCollection(const PositionsVectorCollection& electronCollection);

    explicit ElectronsVectorCollection(const PositionsVectorCollection& electronCollection,
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
