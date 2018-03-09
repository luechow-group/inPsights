//
// Created by Leonard Reuter on 09.03.18.
//

#include "AtomsVectorCollection.h"

AtomsVectorCollection::AtomsVectorCollection()
        : ParticlesVectorCollection(),
          elementTypesVector_(ElementTypesVector()) {}

AtomsVectorCollection::AtomsVectorCollection(const ElementTypesVector &elementTypesVector)
        : ParticlesVectorCollection(),
          elementTypesVector_(elementTypesVector) {}

AtomsVectorCollection::AtomsVectorCollection(const AtomCollection &atomCollection)
        : AtomsVectorCollection(std::vector<AtomCollection>({atomCollection})){}

AtomsVectorCollection::AtomsVectorCollection(const std::vector<AtomCollection> &atomCollectionVector)
        : ParticlesVectorCollection(),
          elementTypesVector_(atomCollectionVector[0].elementTypesVector()) {

    if ( !atomCollectionVector.empty() ){
        for (const auto &atomCollection : atomCollectionVector) {
            positionsVectorCollection_.append(atomCollection.positionsVector());

            assert(elementTypesVector_.elementTypesAsEigenVector()
                   == atomCollection.elementTypesVector().elementTypesAsEigenVector()
                   && "All AtomCollection s must have the same ElementTypesVector.");
        }
    }
}

AtomsVectorCollection::AtomsVectorCollection(const PositionsVectorCollection &positionsVectorCollection)
        : AtomsVectorCollection(positionsVectorCollection,
                              ElementTypesVector(positionsVectorCollection.numberOfPositionsEntities())) {
}

AtomsVectorCollection::AtomsVectorCollection(const PositionsVectorCollection &positionsVectorCollection,
                                         const ElementTypesVector &elementTypesVector)
        : ParticlesVectorCollection(positionsVectorCollection),
          elementTypesVector_(elementTypesVector) {

    assert(numberOfEntities() == positionsVectorCollection_.numberOfEntities()
           && numberOfEntities() == elementTypesVector_.numberOfEntities()
           && "The number of entities in ParticlesVectorCollection, PositionsVectorCollection, and ElementTypesVector must match.");
}

AtomCollection AtomsVectorCollection::operator[](long i) const {
    return AtomCollection(positionsVectorCollection_[i],elementTypesVector_);
}

const ElementTypesVector& AtomsVectorCollection::elementTypesVector() const{
    return elementTypesVector_;
}

ElementTypesVector &AtomsVectorCollection::elementTypesVector() {
    return elementTypesVector_;
}

void AtomsVectorCollection::insert(const AtomCollection &atomCollection, long i) {
    if (elementTypesVector_.numberOfEntities() != 0){
        assert(elementTypesVector_.elementTypesAsEigenVector()
               == atomCollection.elementTypesVector().elementTypesAsEigenVector());
    }
    else{
        elementTypesVector_ = atomCollection.elementTypesVector();
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
        elementTypesVector_.permute(i,j);
    }
}
