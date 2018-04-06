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

ElectronsVectorCollection::ElectronsVectorCollection(const ElectronsVector &electronsVector)
        : ElectronsVectorCollection(std::vector<ElectronsVector>({electronsVector})){}

ElectronsVectorCollection::ElectronsVectorCollection(const std::vector<ElectronsVector> &electronsVectorVector)
        : ParticlesVectorCollection(),
          spinTypesVector_(0) {

    if ( !electronsVectorVector.empty() ){
        spinTypesVector_ = electronsVectorVector[0].spinTypesVector();
        for (const auto &electronsVector : electronsVectorVector) {
            append(electronsVector);
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
           && "The number of entities in ParticlesVectorCollection and PositionsVectorCollection must be equal.");

    assert(positionsVectorCollection.numberOfPositionsEntities() == spinTypesVector_.numberOfEntities()
           && "The number of entities in PositionsVector and SpinTypesVector must be equal.");
}

ElectronsVector ElectronsVectorCollection::operator[](long i) const {
    return ElectronsVector(positionsVectorCollection_[i],spinTypesVector_);
}

const SpinTypesVector& ElectronsVectorCollection::spinTypesVector() const{
    return spinTypesVector_;
}

SpinTypesVector &ElectronsVectorCollection::spinTypesVector() {
    return spinTypesVector_;
}

void ElectronsVectorCollection::insert(const ElectronsVector &electronsVector, long i) {
    if (spinTypesVector_.numberOfEntities() != 0) {
        assert(spinTypesVector_.spinTypesAsEigenVector()
               == electronsVector.spinTypesVector().spinTypesAsEigenVector());
    }
    else{
        spinTypesVector_ = electronsVector.spinTypesVector();
    }
    positionsVectorCollection_.insert(electronsVector.positionsVector(), i);
    incrementNumberOfEntities();
}

void ElectronsVectorCollection::append(const ElectronsVector &electronsVector) {
    insert(electronsVector,numberOfEntities());
}

void ElectronsVectorCollection::prepend(const ElectronsVector &electronsVector) {
    insert(electronsVector,0);
}

void ElectronsVectorCollection::permute(long i, long j) {
    if(i != j) {
        positionsVectorCollection_.permute(i,j);
        spinTypesVector_.permute(i,j);
    }
}
