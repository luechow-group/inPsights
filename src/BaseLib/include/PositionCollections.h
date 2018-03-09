//
// Created by Michael Heuer on 08.03.17.
//

#ifndef AMOLQCGUI_POSITIONCOLLECTIONS_H
#define AMOLQCGUI_POSITIONCOLLECTIONS_H

#include <vector>
#include "PositionCollection.h"

class PositionCollections : public AbstractCollection{
public:
    PositionCollections();
    explicit PositionCollections(const std::vector<PositionCollection> &positionCollections);

    PositionCollection operator[](long i) const;

    void insert (const PositionCollection& positionCollection, long i);
    void append (const PositionCollection& positionCollection);
    void prepend(const PositionCollection& positionCollection);
    void permute(long i, long j) override;


    const std::vector<PositionCollection>& positionCollections() const;
    std::vector<PositionCollection>& positionCollections();

    long numberOfPositionsEntities() const;

private:
    std::vector<PositionCollection> positionCollections_;
    long numberOfPositionEntities_;
};

#endif //AMOLQCGUI_POSITIONCOLLECTIONS_H
