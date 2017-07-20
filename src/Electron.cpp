//
// Created by Morian Sonnet on 16.05.2017.
//

#include "Electron.h"

Electron::Electron(spintype spin, double x, double y, double z): spin(spin), Particle(x,y,z) {
}

void Electron::setAssignedCore(int toAssignCore) {
    this->assignedCore=toAssignCore;
}

spintype Electron::getSpin() const {
    return spin;
}
