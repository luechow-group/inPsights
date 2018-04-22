//
// Created by Michael Heuer on 29.10.17.
//

#include "Electron.h"
#include "ToString.h"

using namespace Eigen;
using namespace Spin;

Electron::Electron(const Eigen::Vector3d& position, const SpinType & spinType)
        : Particle(position, int(spinType)+Spin::storageShift)
{}

Spin::SpinType Electron::spinType() const {
    return static_cast<SpinType>(type_-storageShift);
}

void Electron::setSpinType(const Spin::SpinType & spinType) {
    type_ = int(spinType)+storageShift;
}

int Electron::charge() const {
    return -1;
}

std::ostream& operator<< (std::ostream& os, const Electron& elec) {
    os << elec.toString();
    return os;
}

std::string Electron::toString() const {

    return "e" + Spin::toString(spinType()) + ToString::vector3dToString(position_);
}
