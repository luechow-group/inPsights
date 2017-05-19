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
