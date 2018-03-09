//
// Created by Michael Heuer on 29.10.17.
//

#ifndef AMOLQCGUI_ELECTRON_H
#define AMOLQCGUI_ELECTRON_H

#include "Particle.h"
#include "SpinType.h"

class Electron : public Particle{
public:
    Electron(const Eigen::Vector3d & position, const Spin::SpinType& spinType = Spin::SpinType::none);

    Electron(const Particle& particle, const Spin::SpinType& spinType = Spin::SpinType::none);

    Spin::SpinType spinType()const;

    void setSpinType(const Spin::SpinType & spinType);

    int charge() override;

private:
    Spin::SpinType spinType_;
};

#endif //AMOLQCGUI_ELECTRON_H
