//
// Created by Michael Heuer on 29.10.17.
//

#include "ElectronsVector.h"
#include "ToString.h"

ElectronsVector::ElectronsVector(const Eigen::VectorXd &positions)
        : ParticlesVector(PositionsVector(positions)),
          spinTypesVector_(numberOfEntities()) {
}

ElectronsVector::ElectronsVector(const Eigen::VectorXd &positions, const Eigen::VectorXi &spinTypes)
        : ParticlesVector(PositionsVector(positions)),
          spinTypesVector_(spinTypes) {

    assert(numberOfEntities() == positionsVector_.numberOfEntities()
           && numberOfEntities() == spinTypesVector_.numberOfEntities()
           && "The number of entities in ParticlesVector, PositionsVector, and SpinTypesVector must match.");
}

ElectronsVector::ElectronsVector(const PositionsVector &positionsVector,
                                       const SpinTypesVector &spinTypesVector)
        : ParticlesVector(PositionsVector(positionsVector)),
          spinTypesVector_(spinTypesVector)
{}

Electron ElectronsVector::operator[](long i) const {
    return Electron{positionsVector_[i], spinTypesVector_[i]};
}

void ElectronsVector::insert(const Electron& electron, long i) {
    assert(i >= 0 && "The index must be positive.");
    assert(i <= numberOfEntities() && "The index must be smaller than the number of entities.");

    positionsVector_.insert(electron.position(),i);
    spinTypesVector_.insert(electron.spinType(),i);
    incrementNumberOfEntities();

    assert(positionsVector_.numberOfEntities() == numberOfEntities());
    assert(spinTypesVector_.numberOfEntities() == numberOfEntities());
}

void ElectronsVector::prepend(const Electron& electron) {
    this->insert(electron,0);
}

void ElectronsVector::append(const Electron& electron) {
    this->insert(electron,numberOfEntities());
}

void ElectronsVector::permute(long i, long j) {
    positionsVector_.permute(i,j);
    spinTypesVector_.permute(i,j);
}

const SpinTypesVector &ElectronsVector::spinTypesVector() const {
    return spinTypesVector_;
}

SpinTypesVector &ElectronsVector::spinTypesVector() {
    return spinTypesVector_;
}

std::ostream& operator<<(std::ostream& os, const ElectronsVector& ec){
    for (unsigned long i = 0; i < ec.numberOfEntities(); i++) {
        os << ToString::unsignedLongToString(i + 1) << " " << ec[i] << std::endl;
    }
    return os;
}
