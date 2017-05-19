//
// Created by Moria on 16.05.2017.
//

#include "Core.h"
#include <iostream>
#include "pse.h"

Core::Core(std::string elementType, double x, double y, double z) : Particle(x, y, z), elementType(elementType) {
    std::cout << "The Element Number of the Element " << elementType << " is " << Pse::findElement(elementType) << std::endl;
    std::cout << "The position of this Core is \n" << this->getPosition() << std::endl;
    charge=Pse::findElement(elementType);
}

int Core::getCharge() const {
    return charge;
}

void Core::setCharge(int charge) {
    Core::charge = charge;
}

void Core::setAssignedElectrons(const std::vector<int> &toAssignElectrons) {
    assignedElectrons=toAssignElectrons;
    std::cout << "Assigned the Electrons ";
    for(int i=0;i<assignedElectrons.size();i++){
        std::cout << assignedElectrons[i] << ' ';
    }
    std::cout << std::endl;

}

int Core::getLocalSpinQuantumNumber(const std::vector<Electron> &electrons) {
    int localSpin=0;
    for(std::vector<int>::const_iterator i=assignedElectrons.begin();i!=assignedElectrons.end();i++){
        localSpin+=(int)electrons[*i].getSpin();
    }
    return localSpin;
}
