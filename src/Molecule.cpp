//
// Created by Moria on 19.05.2017.
//

#include "Molecule.h"

Molecule::Molecule() {

}

Molecule::~Molecule() {

}

void Molecule::addCore(std::string &elementType, double x, double y, double z) {
    this->cores.emplace_back(elementType,x,y,z);

}
