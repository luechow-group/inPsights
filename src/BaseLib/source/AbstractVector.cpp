//
// Created by Michael Heuer on 07.03.18.
//

#include "AbstractVector.h"
#include <assert.h>

AbstractVector::AbstractVector(long numberOfEntities)
        : numberOfEntities_(numberOfEntities){
    assert(numberOfEntities >= 0 && "The number of Entities must be positive.");
};

void AbstractVector::incrementNumberOfEntities(){
    numberOfEntities_++;
}

long AbstractVector::numberOfEntities() const{
    return numberOfEntities_;
}

void AbstractVector::setNumberOfEntities(long numberOfEntities){
    assert(numberOfEntities >= 0  && "The number of Entities must be positive.");
    numberOfEntities_ = numberOfEntities;
}

long AbstractVector::calculateIndex(long i) const {
    assert(i < numberOfEntities() && "Index is out of bounds");
    assert(i >= -numberOfEntities() && "Reverse index is out of bounds");
    if (i >= 0) return i;
    return (numberOfEntities()+i);
}
