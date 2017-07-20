//
// Created by Morian Sonnet on 19.05.2017.
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
}

void Molecule::assign(const Assignment &assignment) {
    for(int i=0;i<assignment.size();i++){
        cores[assignment[i].first].setAssignedElectrons(assignment[i].second);
        for(int j=0;j<assignment[i].second.size();j++){
            electrons[assignment[i].second[j]].setAssignedCore(assignment[i].first);
        }
    }
}

int Molecule::getTotalSpinProjectionQuantumNumber() {
    int totalSpin=0;
    for(std::vector<Electron>::const_iterator i=electrons.begin(); i!=electrons.end();i++){
        totalSpin+=(int)(*i).getSpin();
    }
    return totalSpin;
}

int Molecule::getLocalSpinProjectionQuantumNumber(int coreToLookAt) {
    return cores[coreToLookAt].getLocalSpinQuantumNumber(electrons);
}

int Molecule::getLocalSpinProjectionQuantumNumber(std::vector<int> coresToLookAt) {
    if(coresToLookAt.size()==0)return getTotalSpinProjectionQuantumNumber();
    int LocalFragmentSpin=0;
    for(std::vector<int>::const_iterator i=coresToLookAt.begin();i!=coresToLookAt.end();i++){
        LocalFragmentSpin+= cores[*i].getLocalSpinQuantumNumber(electrons);
    }
    return LocalFragmentSpin;
}
