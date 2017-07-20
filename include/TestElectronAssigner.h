//
// Created by Morian Sonneton 19.05.2017.
//

#ifndef LOCALSPINMULTIPLICITY_TESTELECTRONASSIGNER_H
#define LOCALSPINMULTIPLICITY_TESTELECTRONASSIGNER_H
#include "ElectronAssigner.h"

class TestElectronAssigner : public ElectronAssigner{
public:
    virtual Assignment assign(const std::vector<Core> &cores, const std::vector<Particle> &electrons) override;
    virtual Assignment assign(const std::vector<Core> &cores, const std::vector<Electron> &electrons) override;
};


#endif //LOCALSPINMULTIPLICITY_TESTELECTRONASSIGNER_H
