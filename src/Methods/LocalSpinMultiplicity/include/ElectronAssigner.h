//
// Created by Morian Sonnet on 19.05.2017.
//

#ifndef LOCALSPINMULTIPLICITY_ELECTRONASSIGNER_H
#define LOCALSPINMULTIPLICITY_ELECTRONASSIGNER_H
#include "Core.h"
#include "Electron.h"
#include <vector>
#include "Assignment.h"

/*
 * This class represents a ElectronAssigner, assigning Electrons to Cores.
 * This is the abstract base class. Different Implementations are possible.
 * The assign function needs to be able to assign Electrons without Spin (bmax) and Electrons (bstart) to Cores.
 * =>function is doubly implemented. Solution is given by using Particle pointers, since Spin is unused here.
 * This is not implemented yet.
 */
class ElectronAssigner {
public:
    ElectronAssigner();
    virtual ~ElectronAssigner();
    virtual Assignment assign(const std::vector<Core> &, const std::vector<Particle> &) = 0;
    virtual Assignment assign(const std::vector<Core> &, const std::vector<Electron> &) = 0;
private:
};


#endif //LOCALSPINMULTIPLICITY_ELECTRONASSIGNER_H
