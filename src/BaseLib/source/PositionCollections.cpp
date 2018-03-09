//
// Created by Michael Heuer on 30.10.17.
//

#include <ElectronCollection.h>
#include "PositionCollections.h"

PositionCollections::PositionCollections()
        : AbstractCollection(0),
          positionCollections_(0),
          numberOfPositionEntities_(0)
{}

PositionCollections::PositionCollections(const std::vector<PositionCollection> &positionCollections)
        : AbstractCollection(positionCollections.size()),
          numberOfPositionEntities_(0)
{
    if( !positionCollections.empty()) {
        numberOfPositionEntities_ = positionCollections[0].numberOfEntities();

        for (const auto &positionCollection : positionCollections) {
             assert(numberOfPositionEntities_ == positionCollection.numberOfEntities()
                    && "All position collections must contain the same number of positions.");
        }
    }
}

long PositionCollections::numberOfPositionsEntities() const {
    return numberOfPositionEntities_;
}

PositionCollection PositionCollections::operator[](long i) const {
    return positionCollections_[calculateIndex(i)];
}

void PositionCollections::insert(const PositionCollection &positionCollection, long i) {
    assert(i >= 0 && "The index must be positive.");
    assert(i <= numberOfEntities() && "The index must be smaller than the number of entities.");
    if(numberOfEntities() == 0) {
        numberOfPositionEntities_ = positionCollection.numberOfEntities();
    }
    else {
        assert(numberOfPositionEntities_ == positionCollection.numberOfEntities()
               && "All position collections must contain the same number of positions.");
    }

    positionCollections_.insert(positionCollections_.begin()+i,positionCollection);
    incrementNumberOfEntities();
}

void PositionCollections::append(const PositionCollection &positionCollection) {
    insert(positionCollection,numberOfEntities());
}

void PositionCollections::prepend(const PositionCollection &positionCollection) {
    insert(positionCollection,0);
}

void PositionCollections::permute(long i, long j) {
    if(i != j) {
        PositionCollection tempi =positionCollections_[calculateIndex(i)];
        PositionCollection tempj =positionCollections_[calculateIndex(j)];

        positionCollections_.at(calculateIndex(i)) = tempj;
        positionCollections_.at(calculateIndex(j)) = tempi;

        //positionCollections_.( static_cast<unsigned long>(calculateIndex(i)) )
        //        = positionCollections_.at( static_cast<unsigned long>(calculateIndex(j)) );
        //positionCollections_.at( static_cast<unsigned long>(calculateIndex(j)) ) = temp;
    }
}

const std::vector<PositionCollection>& PositionCollections::positionCollections() const {
    return positionCollections_;
}

std::vector<PositionCollection>& PositionCollections::positionCollections() {
    return positionCollections_;
}
