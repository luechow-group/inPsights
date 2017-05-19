//
// Created by Moria on 19.05.2017.
//

#include "Molecule.h"
#include <iostream>

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

void Molecule::assign(ElectronAssigner &electronAssigner) {
    this->assign(electronAssigner.assign(cores,electrons));
}

const std::vector<Core> &Molecule::getCores() const {
    return this->cores;
}

void Molecule::cleanElectrons() {
    electrons.clear();
    std::cout << "Cleared all those filthy Electrons" << std::endl;
}

void Molecule::assign(const Assignation &assignation) {
    for(int i=0;i<assignation.size();i++){
        cores[assignation[i].first].setAssignedElectrons(assignation[i].second);
        for(int j=0;j<assignation[i].second.size();j++){
            electrons[assignation[i].second[j]].setAssignedCore(assignation[i].first);
        }
    }
}
