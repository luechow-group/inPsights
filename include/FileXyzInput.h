//
// Created by Morian Sonnet on 19.05.2017.
//

#ifndef LOCALSPINMULTIPLICITY_FILEXYZINPUT_H
#define LOCALSPINMULTIPLICITY_FILEXYZINPUT_H
#include "InputOutput.h"
#include "Molecule.h"
#include <vector>
#include "Assignment.h"
#include "Core.h"
#include "ElectronAssigner.h"
#include "SpinDeterminer.h"

/*
 * This class is responsible for obtaining the data out of the given .xyz and .ref file.
 * It is derived from the general InputOutput class.
 * When reading Electron Assignment it also calls the functions of given SpinDeterminer and ElectronAssigner to obtain a Electron-Core-Assignment.
 * For bmax-Method it can also read the maximum positions and obtain Electron-Core-Assignments which are saved.
 * These are used when  no ElectronAssigner is given.
 */
class FileXyzInput : public InputOutput{
public:
    FileXyzInput(const std::string &refFilename, const std::string &xyzFilename);
    void readMoleculeCores(Molecule& molecule);
    bool readElectronStructure(Molecule &molecule, const SpinDeterminer &spinDeterminer, ElectronAssigner *ea = nullptr);
    void readElectronCoreAssignments(const std::vector<Core> &cores, ElectronAssigner &ea);
    const std::vector<Assignment> &getAssignments() const;
    void printAssignments();
    virtual ~FileXyzInput();
private:
    std::vector<Assignment> assignments;
};


#endif //LOCALSPINMULTIPLICITY_FILEXYZINPUT_H
