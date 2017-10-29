//
// Created by Michael Heuer on 29.10.17.
//

#ifndef AMOLQCGUI_ELECTRON_H
#define AMOLQCGUI_ELECTRON_H

#include "Particle.h"
#include "SpinType.h"

class Electron : public Particle{
public:
    Electron(const Vector3d & position, Spin::SpinType spinType = Spin::SpinType::none);

    Electron(const Particle& particle, Spin::SpinType spinType = Spin::SpinType::none);

    Spin::SpinType spinType();
    void setSpinType(Spin::SpinType spinType);

private:
    Spin::SpinType spinType_;
};

#endif //AMOLQCGUI_ELECTRON_H
