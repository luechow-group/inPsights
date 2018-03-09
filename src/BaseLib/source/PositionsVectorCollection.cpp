//
// Created by Michael Heuer on 30.10.17.
//

#include <ElectronCollection.h>
#include "PositionsVectorCollection.h"

PositionsVectorCollection::PositionsVectorCollection()
        : AbstractCollection(0),
          positionsVectorCollection_(0),
          numberOfPositionEntities_(0)
{}

PositionsVectorCollection::PositionsVectorCollection(const std::vector<PositionCollection> &positionsVectorCollection)
        : AbstractCollection(positionsVectorCollection.size()),
          numberOfPositionEntities_(0)
{
    if( !positionsVectorCollection.empty()) {
        numberOfPositionEntities_ = positionsVectorCollection[0].numberOfEntities();

        for (const auto &positionCollection : positionsVectorCollection) {
             assert(numberOfPositionEntities_ == positionCollection.numberOfEntities()
                    && "All position collections must contain the same number of positions.");
        }
    }
}

long PositionsVectorCollection::numberOfPositionsEntities() const {
    return numberOfPositionEntities_;
}

PositionCollection PositionsVectorCollection::operator[](long i) const {
    return positionsVectorCollection_[calculateIndex(i)];
}

void PositionsVectorCollection::insert(const PositionCollection &positionCollection, long i) {
    assert(i >= 0 && "The index must be positive.");
    assert(i <= numberOfEntities() && "The index must be smaller than the number of entities.");
    if(numberOfEntities() == 0) {
        numberOfPositionEntities_ = positionCollection.numberOfEntities();
    }
    else {
        assert(numberOfPositionEntities_ == positionCollection.numberOfEntities()
               && "All position collections must contain the same number of positions.");
    }

    positionsVectorCollection_.insert(positionsVectorCollection_.begin()+i,positionCollection);
    incrementNumberOfEntities();
}

void PositionsVectorCollection::append(const PositionCollection &positionCollection) {
    insert(positionCollection,numberOfEntities());
}

void PositionsVectorCollection::prepend(const PositionCollection &positionCollection) {
    insert(positionCollection,0);
}

void PositionsVectorCollection::permute(long i, long j) {
    if(i != j) {
        PositionCollection tempi =positionsVectorCollection_[calculateIndex(i)];
        PositionCollection tempj =positionsVectorCollection_[calculateIndex(j)];

        positionsVectorCollection_.at(calculateIndex(i)) = tempj;
        positionsVectorCollection_.at(calculateIndex(j)) = tempi;

        //positionsVectorCollection_.( static_cast<unsigned long>(calculateIndex(i)) )
        //        = positionsVectorCollection_.at( static_cast<unsigned long>(calculateIndex(j)) );
        //positionsVectorCollection_.at( static_cast<unsigned long>(calculateIndex(j)) ) = temp;
    }
}

const std::vector<PositionCollection>& PositionsVectorCollection::positionsVectorCollection() const {
    return positionsVectorCollection_;
}

std::vector<PositionCollection>& PositionsVectorCollection::positionsVectorCollection() {
    return positionsVectorCollection_;
}
