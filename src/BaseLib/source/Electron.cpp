//
// Created by Michael Heuer on 29.10.17.
//

#include "Electron.h"

Electron::Electron(const Vector3d & position,Spin::SpinType spinType)
        : Particle(position),
          spinType_(spinType)
{}

