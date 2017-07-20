//
// Created by Morian Sonneton 16.05.2017.
//

#ifndef LOCALSPINMULTIPLICITY_ELECTRON_H
#define LOCALSPINMULTIPLICITY_ELECTRON_H
#include "Particle.h"

enum struct spintype {SPIN_ALPHA = 1, SPIN_BETA = -1, SPIN_NA=0};

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
