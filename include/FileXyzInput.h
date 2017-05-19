//
// Created by Moria on 19.05.2017.
//

#ifndef LOCALSPINMULTIPLICITY_FILEXYZINPUT_H
#define LOCALSPINMULTIPLICITY_FILEXYZINPUT_H
#include "InputOutput.h"
#include "Molecule.h"


class FileXyzInput : public InputOutput{
public:
    FileXyzInput(const std::string &refFilename, const std::string &xyzFilename);
    void readMoleculeCores(Molecule& molecule);
    void readElectronStructure(Molecule& molecule);
    virtual ~FileXyzInput();

};


#endif //LOCALSPINMULTIPLICITY_FILEXYZINPUT_H
