//
// Created by Leonard Reuter on 09.03.18.
//

#include "AtomCollections.h"

AtomCollections::AtomCollections()
        : ParticleCollections(),
          elementTypeCollection_(ElementTypeCollection()) {}

AtomCollections::AtomCollections(const ElementTypeCollection &elementTypeCollection)
        : ParticleCollections(),
          elementTypeCollection_(elementTypeCollection) {}

AtomCollections::AtomCollections(const AtomCollection &atomCollection)
        : AtomCollections(std::vector<AtomCollection>({atomCollection})){}

AtomCollections::AtomCollections(const std::vector<AtomCollection> &atomCollectionVector)
        : ParticleCollections(),
          elementTypeCollection_(atomCollectionVector[0].elementTypeCollection()) {

    if ( !atomCollectionVector.empty() ){
        for (const auto &atomCollection : atomCollectionVector) {
            positionCollections_.append(atomCollection.positionCollection());

            assert(elementTypeCollection_.elementTypesAsEigenVector()
                   == atomCollection.elementTypeCollection().elementTypesAsEigenVector()
                   && "All AtomCollection s must have the same ElementTypeCollection.");
        }
    }
}

AtomCollections::AtomCollections(const PositionCollections &positionCollections)
        : AtomCollections(positionCollections,
                              ElementTypeCollection(positionCollections.numberOfPositionsEntities())) {
}

AtomCollections::AtomCollections(const PositionCollections &positionCollections,
                                         const ElementTypeCollection &elementTypeCollection)
        : ParticleCollections(positionCollections),
          elementTypeCollection_(elementTypeCollection) {

    assert(numberOfEntities() == positionCollections_.numberOfEntities()
           && numberOfEntities() == elementTypeCollection_.numberOfEntities()
           && "The number of entities in ParticleCollections, PositionCollections, and ElementTypeCollection must match.");
}

AtomCollection AtomCollections::operator[](long i) const {
    return AtomCollection(positionCollections_[i],elementTypeCollection_);
}

const ElementTypeCollection& AtomCollections::elementTypeCollection() const{
    return elementTypeCollection_;
}

ElementTypeCollection &AtomCollections::elementTypeCollection() {
    return elementTypeCollection_;
}

void AtomCollections::insert(const AtomCollection &atomCollection, long i) {
    if (elementTypeCollection_.numberOfEntities() != 0){
        assert(elementTypeCollection_.elementTypesAsEigenVector()
               == atomCollection.elementTypeCollection().elementTypesAsEigenVector());
    }
    else{
        elementTypeCollection_ = atomCollection.elementTypeCollection();
    }
    positionCollections_.insert(atomCollection.positionCollection(), i);
    incrementNumberOfEntities();
}

void AtomCollections::append(const AtomCollection &atomCollection) {
    insert(atomCollection,numberOfEntities());
}

void AtomCollections::prepend(const AtomCollection &atomCollection) {
    insert(atomCollection,0);
}

void AtomCollections::permute(long i, long j) {
    if(i != j) {
        positionCollections_.permute(i,j);
        elementTypeCollection_.permute(i,j);
    }
}
