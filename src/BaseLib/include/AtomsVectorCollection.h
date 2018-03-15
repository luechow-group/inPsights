//
// Created by Leonard Reuter on 09.03.17.
//

#ifndef AMOLQCPP_ATOMCOLLECTIONS_H
#define AMOLQCPP_ATOMCOLLECTIONS_H

#include "ParticlesVectorCollection.h"
#include "AtomsVector.h"
#include "ElementTypesVector.h"

class AtomsVectorCollection : public ParticlesVectorCollection{
public:
    AtomsVectorCollection();
    explicit AtomsVectorCollection(const ElementTypesVector& elementTypesVector);
    explicit AtomsVectorCollection(const AtomsVector& atomsVector);
    explicit AtomsVectorCollection(const std::vector<AtomsVector>& atomsVectorVector);
    explicit AtomsVectorCollection(const PositionsVectorCollection& atomsVector);

    explicit AtomsVectorCollection(const PositionsVectorCollection& atomsVector,
                                 const ElementTypesVector& elementTypesVector);

    AtomsVector operator[](long i) const;

    const ElementTypesVector& elementTypesVector() const;
    ElementTypesVector& elementTypesVector();

    void insert (const AtomsVector& atomsVector, long i);
    void append (const AtomsVector& atomsVector);
    void prepend(const AtomsVector& atomsVector);
    void permute(long i, long j) override;

private:
    ElementTypesVector elementTypesVector_;
};

#endif //AMOLQCPP_ATOMCOLLECTIONS_H
