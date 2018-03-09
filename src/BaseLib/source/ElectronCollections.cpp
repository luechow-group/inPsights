//
// Created by Michael Heuer on 30.10.17.
//

#include "ElectronCollections.h"

ElectronCollections::ElectronCollections()
        : ParticleCollections(),
          spinTypeCollection_(SpinTypeCollection()) {}

ElectronCollections::ElectronCollections(const SpinTypeCollection &spinTypeCollection)
        : ParticleCollections(),
          spinTypeCollection_(spinTypeCollection) {}

ElectronCollections::ElectronCollections(const ElectronCollection &electronCollection)
        : ElectronCollections(std::vector<ElectronCollection>({electronCollection})){}

ElectronCollections::ElectronCollections(const std::vector<ElectronCollection> &electronCollectionVector)
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

ElectronCollections::ElectronCollections(const PositionCollections &positionCollections)
        : ElectronCollections(positionCollections,
                              SpinTypeCollection(positionCollections.numberOfPositionsEntities())) {
}

ElectronCollections::ElectronCollections(const PositionCollections &positionCollections,
                                         const SpinTypeCollection &spinTypeCollection)
        : ParticleCollections(positionCollections),
          spinTypeCollection_(spinTypeCollection) {

    assert(numberOfEntities() == positionCollections_.numberOfEntities()
           && numberOfEntities() == spinTypeCollection_.numberOfEntities()
           && "The number of entities in ParticleCollections, PositionCollections, and SpinTypeCollection must match.");
}

ElectronCollection ElectronCollections::operator[](long i) const {
    return ElectronCollection(positionCollections_[i],spinTypeCollection_);
}

const SpinTypeCollection& ElectronCollections::spinTypeCollection() const{
    return spinTypeCollection_;
}

SpinTypeCollection &ElectronCollections::spinTypeCollection() {
    return spinTypeCollection_;
}

void ElectronCollections::insert(const ElectronCollection &electronCollection, long i) {
    assert(spinTypeCollection_.spinTypesAsEigenVector()
           == electronCollection.spinTypeCollection().spinTypesAsEigenVector());
    positionCollections_.insert(electronCollection.positionCollection(), i);
    incrementNumberOfEntities();
}

void ElectronCollections::append(const ElectronCollection &electronCollection) {
    insert(electronCollection,numberOfEntities());
}

void ElectronCollections::prepend(const ElectronCollection &electronCollection) {
    insert(electronCollection,0);
}

void ElectronCollections::permute(long i, long j) {
    if(i != j) {
        positionCollections_.permute(i,j);
        spinTypeCollection_.permute(i,j);
    }
}
