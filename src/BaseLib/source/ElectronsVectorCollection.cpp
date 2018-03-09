//
// Created by Michael Heuer on 30.10.17.
//

#include "ElectronsVectorCollection.h"

ElectronsVectorCollection::ElectronsVectorCollection()
        : ParticlesVectorCollection(),
          spinTypesVector_(SpinTypesVector()) {}

ElectronsVectorCollection::ElectronsVectorCollection(const SpinTypesVector &spinTypesVector)
        : ParticlesVectorCollection(),
          spinTypesVector_(spinTypesVector) {}

ElectronsVectorCollection::ElectronsVectorCollection(const ElectronCollection &electronCollection)
        : ElectronsVectorCollection(std::vector<ElectronCollection>({electronCollection})){}

ElectronsVectorCollection::ElectronsVectorCollection(const std::vector<ElectronCollection> &electronCollectionVector)
        : ParticlesVectorCollection(),
          spinTypesVector_(electronCollectionVector[0].spinTypesVector()) {

    if ( !electronCollectionVector.empty() ){
        for (const auto &electronCollection : electronCollectionVector) {
            positionsVectorCollection_.append(electronCollection.positionsVector());

            assert(spinTypesVector_.spinTypesAsEigenVector()
                   == electronCollection.spinTypesVector().spinTypesAsEigenVector()
                   && "All ElectronCollection s must have the same SpinTypesVector.");
        }
    }
}

ElectronsVectorCollection::ElectronsVectorCollection(const PositionsVectorCollection &positionsVectorCollection)
        : ElectronsVectorCollection(positionsVectorCollection,
                              SpinTypesVector(positionsVectorCollection.numberOfPositionsEntities())) {
}

ElectronsVectorCollection::ElectronsVectorCollection(const PositionsVectorCollection &positionsVectorCollection,
                                         const SpinTypesVector &spinTypesVector)
        : ParticlesVectorCollection(positionsVectorCollection),
          spinTypesVector_(spinTypesVector) {

    assert(numberOfEntities() == positionsVectorCollection_.numberOfEntities()
           && numberOfEntities() == spinTypesVector_.numberOfEntities()
           && "The number of entities in ParticlesVectorCollection, PositionsVectorCollection, and SpinTypesVector must match.");
}

ElectronCollection ElectronsVectorCollection::operator[](long i) const {
    return ElectronCollection(positionsVectorCollection_[i],spinTypesVector_);
}

const SpinTypesVector& ElectronsVectorCollection::spinTypesVector() const{
    return spinTypesVector_;
}

SpinTypesVector &ElectronsVectorCollection::spinTypesVector() {
    return spinTypesVector_;
}

void ElectronsVectorCollection::insert(const ElectronCollection &electronCollection, long i) {
    if (spinTypesVector_.numberOfEntities() != 0) {
        assert(spinTypesVector_.spinTypesAsEigenVector()
               == electronCollection.spinTypesVector().spinTypesAsEigenVector());
    }
    else{
        spinTypesVector_ = electronCollection.spinTypesVector();
    }
    positionsVectorCollection_.insert(electronCollection.positionsVector(), i);
    incrementNumberOfEntities();
}

void ElectronsVectorCollection::append(const ElectronCollection &electronCollection) {
    insert(electronCollection,numberOfEntities());
}

void ElectronsVectorCollection::prepend(const ElectronCollection &electronCollection) {
    insert(electronCollection,0);
}

void ElectronsVectorCollection::permute(long i, long j) {
    if(i != j) {
        positionsVectorCollection_.permute(i,j);
        spinTypesVector_.permute(i,j);
    }
}
