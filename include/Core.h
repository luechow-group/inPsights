//
// Created by Moria on 16.05.2017.
//

#ifndef LOCALSPINMULTIPLICITY_CORE_H
#define LOCALSPINMULTIPLICITY_CORE_H
#include "Particle.h"

class Core :public Particle {
public:
    Core(std::string elementType, double x, double y, double z);
private:
    std::string elementType;
 // ULF
};


#endif //LOCALSPINMULTIPLICITY_CORE_H
