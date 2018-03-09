//
// Created by Michael Heuer on 08.03.17.
//

#ifndef AMOLQCGUI_POSITIONCOLLECTIONS_H
#define AMOLQCGUI_POSITIONCOLLECTIONS_H

#include <vector>
#include "PositionCollection.h"

class PositionsVectorCollection : public AbstractCollection{
public:
    PositionsVectorCollection();
    explicit PositionsVectorCollection(const std::vector<PositionCollection> &positionsVectorCollection);

    PositionCollection operator[](long i) const;

    void insert (const PositionCollection& positionCollection, long i);
    void append (const PositionCollection& positionCollection);
    void prepend(const PositionCollection& positionCollection);
    void permute(long i, long j) override;


    const std::vector<PositionCollection>& positionsVectorCollection() const;
    std::vector<PositionCollection>& positionsVectorCollection();

    long numberOfPositionsEntities() const;

private:
    std::vector<PositionCollection> positionsVectorCollection_;
    long numberOfPositionEntities_;
};

#endif //AMOLQCGUI_POSITIONCOLLECTIONS_H
