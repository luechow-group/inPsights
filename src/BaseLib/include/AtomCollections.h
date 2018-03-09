//
// Created by Leonard Reuter on 09.03.17.
//

#ifndef AMOLQCGUI_ATOMCOLLECTIONS_H
#define AMOLQCGUI_ATOMCOLLECTIONS_H

#include "ParticleCollections.h"
#include "AtomCollection.h"
#include "ElementTypeCollection.h"

class AtomCollections : public ParticleCollections{
public:
    AtomCollections();
    explicit AtomCollections(const ElementTypeCollection& elementTypeCollection);
    explicit AtomCollections(const AtomCollection& atomCollection);
    explicit AtomCollections(const std::vector<AtomCollection>& atomCollectionVector);
    explicit AtomCollections(const PositionCollections& atomCollection);

    explicit AtomCollections(const PositionCollections& atomCollection,
                                 const ElementTypeCollection& elementTypeCollection);

    AtomCollection operator[](long i) const;

    const ElementTypeCollection& elementTypeCollection() const;
    ElementTypeCollection& elementTypeCollection();

    void insert (const AtomCollection& atomCollection, long i);
    void append (const AtomCollection& atomCollection);
    void prepend(const AtomCollection& atomCollection);
    void permute(long i, long j) override;

private:
    ElementTypeCollection elementTypeCollection_;
};

#endif //AMOLQCGUI_ATOMCOLLECTIONS_H
