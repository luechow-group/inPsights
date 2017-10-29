//
// Created by Michael Heuer on 29.10.17.
//

#include "Electron.h"

Electron::Electron(const Vector3d & position,Spin::SpinType spinType)
        : Particle(position),
          spinType_(spinType)
{}

Electron::Electron(const Particle &particle, Spin::SpinType spinType)
        : Particle(particle),
          spinType_(spinType)
{}

Spin::SpinType Electron::spin() {
    return spinType_;
}

void Electron::spin(Spin::SpinType spinType) {
    spinType_ = spinType;
}
