//
// Created by Moria on 19.05.2017.
//

#ifndef LOCALSPINMULTIPLICITY_MOLECULE_H
#define LOCALSPINMULTIPLICITY_MOLECULE_H
#include "Core.h"
#include "Electron.h"
#include <vector>
#include "Assignation.h"
#include "ElectronAssigner.h"

class Molecule {
public:
    Molecule();
    virtual ~Molecule();
    void addCore(std::string& elementType,double x, double y, double z);
    void addElectron(spintype spin,double x, double y, double z);
    const std::vector<Core> &getCores() const;
    void cleanElectrons();
    void assign(const Assignation &assignation);
    void assign(ElectronAssigner &electronAssigner);
private:
    std::vector<Core> cores;
    std::vector<Electron> electrons;
};


#endif //LOCALSPINMULTIPLICITY_MOLECULE_H
