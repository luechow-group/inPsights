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

AtomsVectorCollection::AtomsVectorCollection(const AtomsVector &atomsVector)
        : AtomsVectorCollection(std::vector<AtomsVector>({atomsVector})){}

AtomsVectorCollection::AtomsVectorCollection(const std::vector<AtomsVector> &atomsVectorVector)
        : ParticlesVectorCollection(),
          elementTypesVector_(atomsVectorVector[0].elementTypesVector()) {

    if ( !atomsVectorVector.empty() ){
        for (const auto &atomsVector : atomsVectorVector) {
            positionsVectorCollection_.append(atomsVector.positionsVector());

            assert(elementTypesVector_.elementTypesAsEigenVector()
                   == atomsVector.elementTypesVector().elementTypesAsEigenVector()
                   && "All AtomsVector s must have the same ElementTypesVector.");
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

AtomsVector AtomsVectorCollection::operator[](long i) const {
    return AtomsVector(positionsVectorCollection_[i],elementTypesVector_);
}

const ElementTypesVector& AtomsVectorCollection::elementTypesVector() const{
    return elementTypesVector_;
}

ElementTypesVector &AtomsVectorCollection::elementTypesVector() {
    return elementTypesVector_;
}

void AtomsVectorCollection::insert(const AtomsVector &atomsVector, long i) {
    if (elementTypesVector_.numberOfEntities() != 0){
        assert(elementTypesVector_.elementTypesAsEigenVector()
               == atomsVector.elementTypesVector().elementTypesAsEigenVector());
    }
    else{
        elementTypesVector_ = atomsVector.elementTypesVector();
    }
    positionsVectorCollection_.insert(atomsVector.positionsVector(), i);
    incrementNumberOfEntities();
}

void AtomsVectorCollection::append(const AtomsVector &atomsVector) {
    insert(atomsVector,numberOfEntities());
}

void AtomsVectorCollection::prepend(const AtomsVector &atomsVector) {
    insert(atomsVector,0);
}

void AtomsVectorCollection::permute(long i, long j) {
    if(i != j) {
        positionsVectorCollection_.permute(i,j);
        elementTypesVector_.permute(i,j);
    }
}
