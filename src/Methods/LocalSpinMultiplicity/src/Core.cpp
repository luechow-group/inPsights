//
// Created by Morian Sonnet on 16.05.2017.
//

#include "Core.h"
#include <iostream>
#include "pse.h"

Core::Core(std::string elementType, double x, double y, double z) : Particle(x, y, z), elementType(elementType) {
    charge=Pse::findElement(elementType);
}

int Core::getCharge() const {
    return charge;
}

void Core::setCharge(int charge) {      //not used yet, intended for manually changing the amount of assigned electrons
    Core::charge = charge;
}

void Core::setAssignedElectrons(const std::vector<int> &toAssignElectrons) {
    assignedElectrons=toAssignElectrons;
}

int Core::getLocalSpinProjectionQuantumNumber(const std::vector<Electron> &electrons) {
    int localSpin=0;
    for(std::vector<int>::const_iterator i=assignedElectrons.begin();i!=assignedElectrons.end();i++){
        localSpin+=(int)electrons[*i].getSpin();
    }
    return localSpin;
}

const std::string &Core::getElementType() const {
    return elementType;
}
