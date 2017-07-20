//
// Created by Morian Sonneton 16.05.2017.
//

#include "Electron.h"
#include <iostream>

Electron::Electron(spintype spin, double x, double y, double z): spin(spin), Particle(x,y,z) {

    //std::cout << "New Electron created with Spin " << (int)spin << " at position:" << std::endl <<this->getPosition() << std::endl;

}

void Electron::setAssignedCore(int toAssignCore) {
    this->assignedCore=toAssignCore;
    //std::cout<<"Assigned Electron to Core " << assignedCore << std::endl;
}

spintype Electron::getSpin() const {
    return spin;
}
