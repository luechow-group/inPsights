//
// Created by Leonard Reuter on 09.03.18.
//

#include "AtomsVectorCollection.h"

AtomsVectorCollection::AtomsVectorCollection()
        : ParticleCollections(),
          elementTypeCollection_(ElementTypeCollection()) {}

AtomsVectorCollection::AtomsVectorCollection(const ElementTypeCollection &elementTypeCollection)
        : ParticleCollections(),
          elementTypeCollection_(elementTypeCollection) {}

AtomsVectorCollection::AtomsVectorCollection(const AtomCollection &atomCollection)
        : AtomsVectorCollection(std::vector<AtomCollection>({atomCollection})){}

AtomsVectorCollection::AtomsVectorCollection(const std::vector<AtomCollection> &atomCollectionVector)
        : ParticleCollections(),
          elementTypeCollection_(atomCollectionVector[0].elementTypeCollection()) {

    if ( !atomCollectionVector.empty() ){
        for (const auto &atomCollection : atomCollectionVector) {
            positionsVectorCollection_.append(atomCollection.positionsVector());

            assert(elementTypeCollection_.elementTypesAsEigenVector()
                   == atomCollection.elementTypeCollection().elementTypesAsEigenVector()
                   && "All AtomCollection s must have the same ElementTypeCollection.");
        }
    }
}

AtomsVectorCollection::AtomsVectorCollection(const PositionsVectorCollection &positionsVectorCollection)
        : AtomsVectorCollection(positionsVectorCollection,
                              ElementTypeCollection(positionsVectorCollection.numberOfPositionsEntities())) {
}

AtomsVectorCollection::AtomsVectorCollection(const PositionsVectorCollection &positionsVectorCollection,
                                         const ElementTypeCollection &elementTypeCollection)
        : ParticleCollections(positionsVectorCollection),
          elementTypeCollection_(elementTypeCollection) {

    assert(numberOfEntities() == positionsVectorCollection_.numberOfEntities()
           && numberOfEntities() == elementTypeCollection_.numberOfEntities()
           && "The number of entities in ParticleCollections, PositionsVectorCollection, and ElementTypeCollection must match.");
}

AtomCollection AtomsVectorCollection::operator[](long i) const {
    return AtomCollection(positionsVectorCollection_[i],elementTypeCollection_);
}

const ElementTypeCollection& AtomsVectorCollection::elementTypeCollection() const{
    return elementTypeCollection_;
}

ElementTypeCollection &AtomsVectorCollection::elementTypeCollection() {
    return elementTypeCollection_;
}

void AtomsVectorCollection::insert(const AtomCollection &atomCollection, long i) {
    if (elementTypeCollection_.numberOfEntities() != 0){
        assert(elementTypeCollection_.elementTypesAsEigenVector()
               == atomCollection.elementTypeCollection().elementTypesAsEigenVector());
    }
    else{
        elementTypeCollection_ = atomCollection.elementTypeCollection();
    }
    positionsVectorCollection_.insert(atomCollection.positionsVector(), i);
    incrementNumberOfEntities();
}

void AtomsVectorCollection::append(const AtomCollection &atomCollection) {
    insert(atomCollection,numberOfEntities());
}

void AtomsVectorCollection::prepend(const AtomCollection &atomCollection) {
    insert(atomCollection,0);
}

void AtomsVectorCollection::permute(long i, long j) {
    if(i != j) {
        positionsVectorCollection_.permute(i,j);
        elementTypeCollection_.permute(i,j);
    }
}
