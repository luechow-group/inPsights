//
// Created by Michael Heuer on 30.10.17.
//

#include "ElectronsVectorCollection.h"

ElectronsVectorCollection::ElectronsVectorCollection()
        : ParticleCollections(),
          spinTypeCollection_(SpinTypeCollection()) {}

ElectronsVectorCollection::ElectronsVectorCollection(const SpinTypeCollection &spinTypeCollection)
        : ParticleCollections(),
          spinTypeCollection_(spinTypeCollection) {}

ElectronsVectorCollection::ElectronsVectorCollection(const ElectronCollection &electronCollection)
        : ElectronsVectorCollection(std::vector<ElectronCollection>({electronCollection})){}

ElectronsVectorCollection::ElectronsVectorCollection(const std::vector<ElectronCollection> &electronCollectionVector)
        : ParticleCollections(),
          spinTypeCollection_(electronCollectionVector[0].spinTypeCollection()) {

    if ( !electronCollectionVector.empty() ){
        for (const auto &electronCollection : electronCollectionVector) {
            positionCollections_.append(electronCollection.positionCollection());

            assert(spinTypeCollection_.spinTypesAsEigenVector()
                   == electronCollection.spinTypeCollection().spinTypesAsEigenVector()
                   && "All ElectronCollection s must have the same SpinTypeCollection.");
        }
    }
}

ElectronsVectorCollection::ElectronsVectorCollection(const PositionCollections &positionCollections)
        : ElectronsVectorCollection(positionCollections,
                              SpinTypeCollection(positionCollections.numberOfPositionsEntities())) {
}

ElectronsVectorCollection::ElectronsVectorCollection(const PositionCollections &positionCollections,
                                         const SpinTypeCollection &spinTypeCollection)
        : ParticleCollections(positionCollections),
          spinTypeCollection_(spinTypeCollection) {

    assert(numberOfEntities() == positionCollections_.numberOfEntities()
           && numberOfEntities() == spinTypeCollection_.numberOfEntities()
           && "The number of entities in ParticleCollections, PositionCollections, and SpinTypeCollection must match.");
}

ElectronCollection ElectronsVectorCollection::operator[](long i) const {
    return ElectronCollection(positionCollections_[i],spinTypeCollection_);
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
    positionCollections_.insert(electronCollection.positionCollection(), i);
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
        positionCollections_.permute(i,j);
        spinTypeCollection_.permute(i,j);
    }
}
