//
// Created by Moria on 19.05.2017.
//

#ifndef LOCALSPINMULTIPLICITY_MOLECULE_H
#define LOCALSPINMULTIPLICITY_MOLECULE_H
#include "Core.h"
#include <vector>

class Molecule {
public:
    Molecule();
    virtual ~Molecule();
    void addCore(std::string& elementType,double x, double y, double z);
private:
    std::vector<Core> cores;
};


#endif //LOCALSPINMULTIPLICITY_MOLECULE_H
