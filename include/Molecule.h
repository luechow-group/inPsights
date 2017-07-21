//
// Created by Morian Sonnet on 19.05.2017.
//

#ifndef LOCALSPINMULTIPLICITY_MOLECULE_H
#define LOCALSPINMULTIPLICITY_MOLECULE_H
#include "Core.h"
#include "Electron.h"
#include <vector>
#include "Assignment.h"
#include "ElectronAssigner.h"

/*
 * This class represent a Molecule.
 * It contains a set of Cores and Electrons.
 * It can assign its electrons to cores based on a given assignment or given ElectronAssigner.
 * Furthermore it can return the SpinProjectionQuantumNumber of a Fragment, the whole molecule, or a single Atom.
 */
class Molecule {
public:
    Molecule();
    virtual ~Molecule();
    void addCore(std::string& elementType,double x, double y, double z);
    void addElectron(spintype spin,double x, double y, double z);
    const std::vector<Core> &getCores() const;
    void cleanElectrons();
    void assign(const Assignment &assignation);
    void assign(ElectronAssigner &electronAssigner);
    int getTotalSpinProjectionQuantumNumber();
    int getLocalSpinProjectionQuantumNumber(int coreToLookAt);
    int getLocalSpinProjectionQuantumNumber(std::vector<int> coresToLookAt);
private:
    std::vector<Core> cores;
    std::vector<Electron> electrons;
};


#endif //LOCALSPINMULTIPLICITY_MOLECULE_H
