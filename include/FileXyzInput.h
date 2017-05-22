//
// Created by Moria on 19.05.2017.
//

#ifndef LOCALSPINMULTIPLICITY_FILEXYZINPUT_H
#define LOCALSPINMULTIPLICITY_FILEXYZINPUT_H
#include "InputOutput.h"
#include "Molecule.h"
#include <vector>
#include "Assignation.h"
#include "Core.h"
#include "ElectronAssigner.h"
#include "SpinDeterminer.h"


class FileXyzInput : public InputOutput{
public:
    FileXyzInput(const std::string &refFilename, const std::string &xyzFilename);
    void readMoleculeCores(Molecule& molecule);
    int readElectronStructure(Molecule &molecule, const SpinDeterminer &spinDeterminer, ElectronAssigner *ea=0);
    void readElectronCoreAssignations(const std::vector<Core> &cores, ElectronAssigner &ea);
    const std::vector<Assignation> &getAssignations() const;
    void printAssignations();
    virtual ~FileXyzInput();
private:
    std::vector<Assignation> assignations;
};


#endif //LOCALSPINMULTIPLICITY_FILEXYZINPUT_H
