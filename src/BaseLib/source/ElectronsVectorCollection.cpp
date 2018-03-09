//
// Created by Michael Heuer on 30.10.17.
//

#include "ElectronsVectorCollection.h"

ElectronsVectorCollection::ElectronsVectorCollection()
        : ParticlesVectorCollection(),
          spinTypeCollection_(SpinTypeCollection()) {}

ElectronsVectorCollection::ElectronsVectorCollection(const SpinTypeCollection &spinTypeCollection)
        : ParticlesVectorCollection(),
          spinTypeCollection_(spinTypeCollection) {}

ElectronsVectorCollection::ElectronsVectorCollection(const ElectronCollection &electronCollection)
        : ElectronsVectorCollection(std::vector<ElectronCollection>({electronCollection})){}

ElectronsVectorCollection::ElectronsVectorCollection(const std::vector<ElectronCollection> &electronCollectionVector)
        : ParticlesVectorCollection(),
          spinTypeCollection_(electronCollectionVector[0].spinTypeCollection()) {

    if ( !electronCollectionVector.empty() ){
        for (const auto &electronCollection : electronCollectionVector) {
            positionsVectorCollection_.append(electronCollection.positionsVector());

            assert(spinTypeCollection_.spinTypesAsEigenVector()
                   == electronCollection.spinTypeCollection().spinTypesAsEigenVector()
                   && "All ElectronCollection s must have the same SpinTypeCollection.");
        }
    }
}

ElectronsVectorCollection::ElectronsVectorCollection(const PositionsVectorCollection &positionsVectorCollection)
        : ElectronsVectorCollection(positionsVectorCollection,
                              SpinTypeCollection(positionsVectorCollection.numberOfPositionsEntities())) {
}

ElectronsVectorCollection::ElectronsVectorCollection(const PositionsVectorCollection &positionsVectorCollection,
                                         const SpinTypeCollection &spinTypeCollection)
        : ParticlesVectorCollection(positionsVectorCollection),
          spinTypeCollection_(spinTypeCollection) {

    assert(numberOfEntities() == positionsVectorCollection_.numberOfEntities()
           && numberOfEntities() == spinTypeCollection_.numberOfEntities()
           && "The number of entities in ParticlesVectorCollection, PositionsVectorCollection, and SpinTypeCollection must match.");
}

ElectronCollection ElectronsVectorCollection::operator[](long i) const {
    return ElectronCollection(positionsVectorCollection_[i],spinTypeCollection_);
}

const SpinTypeCollection& ElectronsVectorCollection::spinTypeCollection() const{
    return spinTypeCollection_;
}

SpinTypeCollection &ElectronsVectorCollection::spinTypeCollection() {
    return spinTypeCollection_;
}

void ElectronsVectorCollection::insert(const ElectronCollection &electronCollection, long i) {
    if (spinTypeCollection_.numberOfEntities() != 0) {
        assert(spinTypeCollection_.spinTypesAsEigenVector()
               == electronCollection.spinTypeCollection().spinTypesAsEigenVector());
    }
    else{
        spinTypeCollection_ = electronCollection.spinTypeCollection();
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
        spinTypeCollection_.permute(i,j);
    }
}
