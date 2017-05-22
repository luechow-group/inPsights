//
// Created by Moria on 16.05.2017.
//

#ifndef LOCALSPINMULTIPLICITY_CORE_H
#define LOCALSPINMULTIPLICITY_CORE_H

#include <vector>
#include "Particle.h"
#include "Electron.h"

class Core :public Particle {
public:
    Core(std::string elementType, double x, double y, double z);
    void setAssignedElectrons(const std::vector<int> & toAssignElectrons);
private:
    std::string elementType;
    int charge;
public:
    int getCharge() const;
    const std::string &getElementType() const;
    void setCharge(int charge);
    std::vector<int> assignedElectrons;
    int getLocalSpinQuantumNumber(const std::vector<Electron> &electrons);
    // ULF
};


#endif //LOCALSPINMULTIPLICITY_CORE_H
