//
// Created by Moria on 19.05.2017.
//

#include "FileXyzInput.h"
#include "pse.h"
#include "Molecule.h"
#include <iostream>
#include <limits>

FileXyzInput::FileXyzInput(const std::string &refFilename, const std::string &xyzFilename) : InputOutput(refFilename, true) {
    this->openFile(xyzFilename, true);
}

void FileXyzInput::readMoleculeCores(Molecule &molecule) {
    int numCores;
    this->streams[0]>>numCores;
    std::cout << "Num of Cores is " << numCores << std::endl;
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

void FileXyzInput::readElectronStructure(Molecule &molecule) {

}
