//
// Created by Michael Heuer on 29.10.17.
//

#include "Electron.h"
#include "ToString.h"

using namespace Eigen;
using namespace Spin;

Electron::Electron(const Vector3d & position, const SpinType & spinType)
        : Particle(position),
          spinType_(spinType)
{}

Electron::Electron(const Particle &particle, const SpinType & spinType)
        : Particle(particle),
          spinType_(spinType)
{}

Spin::SpinType Electron::spinType() const {
    return spinType_;
}

void Electron::setSpinType(const Spin::SpinType & spinType) {
    spinType_ = spinType;
}

int Electron::charge() const {
    return -1;
}

std::ostream& operator<< (std::ostream& os, const Electron& elec) {
    os << elec.toString();
    return os;
}

std::string Electron::toString() const {

    return "e" + Spin::toString(spinType_) + ToString::vector3dToString(position_);
}
