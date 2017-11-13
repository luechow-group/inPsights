//
// Created by Morian Sonnet on 16.05.2017.
//

#ifndef LOCALSPINMULTIPLICITY_CORE_H
#define LOCALSPINMULTIPLICITY_CORE_H

#include <vector>
#include "Particle.h"
#include "Electron.h"

/*
 * This class represents a Core, being derived from a Particle.
 * Additional features are ...
 * ... a set of assigned Electrons.
 * ... charge of the Core itself.
 * ... element Type.
 */
class Core:public Particle {
private:
    std::string elementType;
    int charge;
    std::vector<int> assignedElectrons;
public:
    Core(std::string elementType, double x, double y, double z);
    int getCharge() const;
    const std::string &getElementType() const;
    void setCharge(int charge);
    void setAssignedElectrons(const std::vector<int> & toAssignElectrons);
    int getLocalSpinProjectionQuantumNumber(const std::vector<Electron> &electrons);
};


#endif //LOCALSPINMULTIPLICITY_CORE_H
