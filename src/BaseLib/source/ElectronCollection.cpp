//
// Created by Michael Heuer on 29.10.17.
//

#include "ElectronCollection.h"
#include "ToString.h"

ElectronCollection::ElectronCollection(const Eigen::VectorXd &positions)
        : ParticlesVector(PositionsVector(positions)),
          spinTypesVector_(numberOfEntities()) {
}

ElectronCollection::ElectronCollection(const Eigen::VectorXd &positions, const Eigen::VectorXi &spinTypes)
        : ParticlesVector(PositionsVector(positions)),
          spinTypesVector_(spinTypes) {

    assert(numberOfEntities() == positionsVector_.numberOfEntities()
           && numberOfEntities() == spinTypesVector_.numberOfEntities()
           && "The number of entities in ParticlesVector, PositionsVector, and SpinTypesVector must match.");
}

ElectronCollection::ElectronCollection(const PositionsVector &positionsVector,
                                       const SpinTypesVector &spinTypesVector)
        : ParticlesVector(PositionsVector(positionsVector)),
          spinTypesVector_(spinTypesVector)
{}

Electron ElectronCollection::operator[](long i) const {
    return Electron{positionsVector_[i], spinTypesVector_[i]};
}

void ElectronCollection::insert(const Electron& electron, long i) {
    assert(i >= 0 && "The index must be positive.");
    assert(i <= numberOfEntities() && "The index must be smaller than the number of entities.");

    positionsVector_.insert(electron.position(),i);
    spinTypesVector_.insert(electron.spinType(),i);
    incrementNumberOfEntities();

    assert(positionsVector_.numberOfEntities() == numberOfEntities());
    assert(spinTypesVector_.numberOfEntities() == numberOfEntities());
}

void ElectronCollection::prepend(const Electron& electron) {
    this->insert(electron,0);
}

void ElectronCollection::append(const Electron& electron) {
    this->insert(electron,numberOfEntities());
}

void ElectronCollection::permute(long i, long j) {
    positionsVector_.permute(i,j);
    spinTypesVector_.permute(i,j);
}

const SpinTypesVector &ElectronCollection::spinTypesVector() const {
    return spinTypesVector_;
}

SpinTypesVector &ElectronCollection::spinTypesVector() {
    return spinTypesVector_;
}

std::ostream& operator<<(std::ostream& os, const ElectronCollection& ec){
    for (unsigned long i = 0; i < ec.numberOfEntities(); i++) {
        os << ToString::intToString(i + 1) << " " << ec[i] << std::endl;
    }
    return os;
}
