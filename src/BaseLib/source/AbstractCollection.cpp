//
// Created by Michael Heuer on 07.03.18.
//

#include "AbstractCollection.h"
#include <assert.h>

AbstractCollection::AbstractCollection(long numberOfEntities)
        : numberOfEntities_(numberOfEntities){
    assert(numberOfEntities >= 0 && "The number of Entities must be positive.");
};

void AbstractCollection::incrementNumberOfEntities(){
    numberOfEntities_++;
}

long AbstractCollection::numberOfEntities() const{
    return numberOfEntities_;
}

void AbstractCollection::setNumberOfEntities(long numberOfEntities){
    assert(numberOfEntities >= 0  && "The number of Entities must be positive.");
    numberOfEntities_ = numberOfEntities;
}

long AbstractCollection::calculateIndex(long i) const {
    assert(i < numberOfEntities() && "Index is out of bounds");
    assert(i >= -numberOfEntities() && "Reverse index is out of bounds");
    if (i >= 0) return i;
    return (numberOfEntities()+i);
}
