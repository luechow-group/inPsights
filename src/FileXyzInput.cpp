//
// Created by Morian Sonnet on 19.05.2017.
//

#include "FileXyzInput.h"
#include "Molecule.h"
#include <iostream>
#include <SpinDeterminer.h>

FileXyzInput::FileXyzInput(const std::string &refFilename, const std::string &xyzFilename) : InputOutput(refFilename, true) {
    this->openFile(xyzFilename, true);
}

void FileXyzInput::readMoleculeCores(Molecule &molecule) {
    int numCores;
    this->streams[0]>>numCores;
    streams[0].ignore(std::numeric_limits<std::streamsize>::max(),'\n');  // go to next line
    for(int i=0;i<numCores;i++){
        int id;
        std::string elementType;
        double x,y,z;
        streams[0]>>id>>elementType>>x>>y>>z;
        molecule.addCore(elementType,x,y,z);
    }


}

FileXyzInput::~FileXyzInput() {
    this->closeAllFiles();

}

int FileXyzInput::readElectronStructure(Molecule &molecule, const SpinDeterminer &spinDeterminer, ElectronAssigner *ea) {
    molecule.cleanElectrons();
    std::string helpString;
    streams[1]>>helpString;
    while(helpString.compare("xyz:")){streams[1]>>helpString;if(streams[1].eof())return 1;}
    int assignmentToUse;
    streams[1]>>assignmentToUse;
    assignmentToUse--;
    streams[1].ignore(std::numeric_limits<std::streamsize>::max(),'\n');  // go to next line
    int numElectrons;
    streams[1]>>numElectrons;
    double x,y,z;
    int index;
    for(int i=0;i<numElectrons;i++){
        streams[1]>>x>>y>>z>>index;
        molecule.addElectron(spinDeterminer.determineSpin(index),x,y,z);
    }
    if(ea){
        molecule.assign(*ea);
    } else {
        molecule.assign(assignments[assignmentToUse]);
    }
    return 0;
}

void FileXyzInput::readElectronCoreAssignments(const std::vector<Core> &cores, ElectronAssigner &ea) {
    int numRefs;
    streams[0]>>numRefs;
    streams[0].ignore(std::numeric_limits<std::streamsize>::max(),'\n');  // go to next line
    for(int i=0;i<numRefs;i++){
        streams[0].ignore(std::numeric_limits<std::streamsize>::max(),'\n');  // go to next line
        std::vector<Particle> tempElectrons;
        int numElectrons;
        streams[0]>>numElectrons;
        for(int j=0;j<numElectrons;j++){
            double x,y,z;
            streams[0]>>x>>y>>z;
            tempElectrons.emplace_back(x,y,z);
        }
        assignments.push_back(ea.assign(cores,tempElectrons));
        tempElectrons.clear();
        streams[0].ignore(std::numeric_limits<std::streamsize>::max(),'\n');  // go to next line
    }
}

void FileXyzInput::printAssignments() {
    for(int i=0;i<assignments.size();i++){
        std::cout << "Printing Assignment number " << i << std::endl;
        for(int j=0;j<assignments[i].size();j++){
            std::cout << "For core " << assignments[i][j].first << " the following electrons are assigned" << std::endl;
            for(int k=0;k<assignments[i][j].second.size();k++)std::cout << assignments[i][j].second[k] << " ";
            std::cout << std::endl;
        }
    }
}

const std::vector<Assignment> &FileXyzInput::getAssignments() const {
    return assignments;
}
