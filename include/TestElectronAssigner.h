//
// Created by Moria on 19.05.2017.
//

#ifndef LOCALSPINMULTIPLICITY_TESTELECTRONASSIGNER_H
#define LOCALSPINMULTIPLICITY_TESTELECTRONASSIGNER_H
#include "ElectronAssigner.h"

class TestElectronAssigner : public ElectronAssigner{
public:
    virtual Assignation assign(const std::vector<Core> &cores, const std::vector<Particle> &electrons) override;
};


#endif //LOCALSPINMULTIPLICITY_TESTELECTRONASSIGNER_H
