//
// Created by Michael Heuer on 29.10.17.
//

#include "ElectronCollection.h"
#include "ToString.h"

ElectronCollection::ElectronCollection(const Eigen::VectorXd &positions)
        : ParticleCollection(PositionCollection(positions)),
          spinTypeCollection_(numberOfEntities()) {
}

ElectronCollection::ElectronCollection(const Eigen::VectorXd &positions, const Eigen::VectorXi &spinTypes)
        : ParticleCollection(PositionCollection(positions)),
          spinTypeCollection_(spinTypes) {

    assert(numberOfEntities() == positionCollection_.numberOfEntities()
           && numberOfEntities() == spinTypeCollection_.numberOfEntities()
           && "The number of entities in ParticleCollection, PositionCollection, and SpinTypeCollection must match.");
}

ElectronCollection::ElectronCollection(const PositionCollection &positionCollection,
                                       const SpinTypeCollection &spinTypeCollection)
        : ParticleCollection(PositionCollection(positionCollection)),
          spinTypeCollection_(spinTypeCollection)
{}

Electron ElectronCollection::operator[](long i) const {
    return Electron{positionCollection_[i], spinTypeCollection_[i]};
}

void ElectronCollection::insert(const Electron& electron, long i) {
    assert(i >= 0 && "The index must be positive.");
    assert(i <= numberOfEntities() && "The index must be smaller than the number of entities.");

    positionCollection_.insert(electron.position(),i);
    spinTypeCollection_.insert(electron.spinType(),i);
    incrementNumberOfEntities();

    assert(positionCollection_.numberOfEntities() == numberOfEntities());
    assert(spinTypeCollection_.numberOfEntities() == numberOfEntities());
}

void ElectronCollection::prepend(const Electron& electron) {
    this->insert(electron,0);
}

void ElectronCollection::append(const Electron& electron) {
    this->insert(electron,numberOfEntities());
}

void ElectronCollection::permute(long i, long j) {
    positionCollection_.permute(i,j);
    spinTypeCollection_.permute(i,j);
}

const SpinTypeCollection &ElectronCollection::spinTypeCollection() const {
    return spinTypeCollection_;
}

SpinTypeCollection &ElectronCollection::spinTypeCollection() {
    return spinTypeCollection_;
}

std::ostream& operator<<(std::ostream& os, const ElectronCollection& ec){
    for (unsigned long i = 0; i < ec.numberOfEntities(); i++) {
        os << ToString::int2string(i+1) << " " << ec[i] << std::endl;
    }
    return os;
}
