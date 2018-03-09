//
// Created by Leonard Reuter on 09.03.17.
//

#ifndef AMOLQCGUI_ATOMCOLLECTIONS_H
#define AMOLQCGUI_ATOMCOLLECTIONS_H

#include "ParticlesVectorCollection.h"
#include "AtomCollection.h"
#include "ElementTypesVector.h"

class AtomsVectorCollection : public ParticlesVectorCollection{
public:
    AtomsVectorCollection();
    explicit AtomsVectorCollection(const ElementTypesVector& elementTypesVector);
    explicit AtomsVectorCollection(const AtomCollection& atomCollection);
    explicit AtomsVectorCollection(const std::vector<AtomCollection>& atomCollectionVector);
    explicit AtomsVectorCollection(const PositionsVectorCollection& atomCollection);

    explicit AtomsVectorCollection(const PositionsVectorCollection& atomCollection,
                                 const ElementTypesVector& elementTypesVector);

    AtomCollection operator[](long i) const;

    const ElementTypesVector& elementTypesVector() const;
    ElementTypesVector& elementTypesVector();

    void insert (const AtomCollection& atomCollection, long i);
    void append (const AtomCollection& atomCollection);
    void prepend(const AtomCollection& atomCollection);
    void permute(long i, long j) override;

private:
    ElementTypesVector elementTypesVector_;
};

#endif //AMOLQCGUI_ATOMCOLLECTIONS_H
