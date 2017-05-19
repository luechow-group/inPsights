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

void Molecule::addElectron(spintype spin, double x, double y, double z) {
    this->electrons.emplace_back(spin,x,y,z);
}

const std::vector<Core> &Molecule::getCores() const {
    return this->cores;
}
