//
// Created by Moria on 16.05.2017.
//

#ifndef LOCALSPINMULTIPLICITY_ELECTRON_H
#define LOCALSPINMULTIPLICITY_ELECTRON_H
#include "Particle.h"

enum struct spintype {SPIN_ALPHA = 1, SPIN_BETA = -1, SPIN_NA=0};

class Electron : public Particle {
public:
    Electron(spintype spin, double x, double y, double z);
private:
    spintype spin;

};


#endif //LOCALSPINMULTIPLICITY_ELECTRON_H
