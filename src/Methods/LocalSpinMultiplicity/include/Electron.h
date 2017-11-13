//
// Created by Morian Sonnet on 16.05.2017.
//

#ifndef LOCALSPINMULTIPLICITY_ELECTRON_H
#define LOCALSPINMULTIPLICITY_ELECTRON_H
#include "Particle.h"

enum struct spintype {SPIN_ALPHA = 1, SPIN_BETA = -1, SPIN_NA=0};

/*
 * This class represents an Electron, derived from a Particle.
 * Additional features are ...
 * ... Spin, which can be alpha or beta (or NA, which is not used yet).
 * ... core, which the electron is assigned to.
 */
class Electron : public Particle {
public:
    Electron(spintype spin, double x, double y, double z);
    void setAssignedCore(int toAssignCore);
    spintype getSpin() const;
private:
    spintype spin;
    int assignedCore=-1;

};


#endif //LOCALSPINMULTIPLICITY_ELECTRON_H
