//
// Created by Michael Heuer on 08.03.17.
//

#ifndef AMOLQCPP_POSITIONCOLLECTIONS_H
#define AMOLQCPP_POSITIONCOLLECTIONS_H

#include <vector>
#include "PositionsVector.h"

class PositionsVectorCollection : public AbstractVector{
public:
    PositionsVectorCollection();
    explicit PositionsVectorCollection(const std::vector<PositionsVector> &positionsVectorCollection);

    PositionsVector operator[](long i) const;

    void insert (const PositionsVector& positionsVector, long i);
    void append (const PositionsVector& positionsVector);
    void prepend(const PositionsVector& positionsVector);
    void permute(long i, long j) override;


    const std::vector<PositionsVector>& positionsVectorCollection() const;
    std::vector<PositionsVector>& positionsVectorCollection();

    long numberOfPositionsEntities() const;

private:
    std::vector<PositionsVector> positionsVectorCollection_;
    long numberOfPositionEntities_;
};

#endif //AMOLQCPP_POSITIONCOLLECTIONS_H
