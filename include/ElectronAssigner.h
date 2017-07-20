//
// Created by Morian Sonnet on 19.05.2017.
//

#ifndef LOCALSPINMULTIPLICITY_ELECTRONASSIGNER_H
#define LOCALSPINMULTIPLICITY_ELECTRONASSIGNER_H
#include "Core.h"
#include "Electron.h"
#include <vector>
#include "Assignment.h"

class ElectronAssigner {
public:
    ElectronAssigner();
    virtual Assignment assign(const std::vector<Core> &, const std::vector<Particle> &) = 0;
    virtual Assignment assign(const std::vector<Core> &, const std::vector<Electron> &) = 0;
    virtual ~ElectronAssigner();
private:
};


#endif //LOCALSPINMULTIPLICITY_ELECTRONASSIGNER_H
