//
// Created by Michael Heuer on 29.10.17.
//

#include "Electron.h"

using namespace Eigen;
using namespace Spin;

Electron::Electron(const Vector3d & position,SpinType spinType)
        : Particle(position),
          spinType_(spinType)
{}

Electron::Electron(const Particle &particle, SpinType spinType)
        : Particle(particle),
          spinType_(spinType)
{}

Spin::SpinType Electron::spinType() const {
    return spinType_;
}
