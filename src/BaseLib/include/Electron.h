//
// Created by Michael Heuer on 29.10.17.
//

#ifndef AMOLQCGUI_ELECTRON_H
#define AMOLQCGUI_ELECTRON_H

#include "Particle.h"
#include "SpinType.h"

class Electron : Particle{
public:
    Electron(const Vector3d & position, Spin::SpinType spinType = Spin::SpinType::None);

private:
    Spin::SpinType spinType_;
};

#endif //AMOLQCGUI_ELECTRON_H
