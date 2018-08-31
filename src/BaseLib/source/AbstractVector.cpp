//
// Created by Michael Heuer on 07.03.18.
//

#include <AbstractVector.h>
#include <assert.h>

AbstractVector::AbstractVector(long numberOfEntities, long entityLength)
        :
        entityLength_(entityLength),
        numberOfEntities_(numberOfEntities){
    assert(numberOfEntities >= 0 && "The number of Entities must be non-negative.");
};

void AbstractVector::incrementNumberOfEntities(){
    numberOfEntities_++;
}

long AbstractVector::numberOfEntities() const{
    return numberOfEntities_;
}

void AbstractVector::setNumberOfEntities(long numberOfEntities){
    assert(numberOfEntities >= 0  && "The number of Entities must be non-negative.");
    numberOfEntities_ = numberOfEntities;
}

long AbstractVector::calculateIndex(long i) const {
    assert(i <= numberOfEntities() && "Index is out of bounds");// less or equal because of insert
    assert(i >= -numberOfEntities() && "Reverse index is out of bounds");
    if (i >= 0) return i*entityLength_;
    return (numberOfEntities()+i)*entityLength_;
}
